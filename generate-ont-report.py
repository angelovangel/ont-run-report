#!/usr/bin/env python3
"""
generate_report.py
==================
Generate an interactive HTML report for ONT PromethION sequencing runs
from a set of pore_activity and throughput CSV files.

Usage
-----
  python generate_report.py [OPTIONS]

Options
-------
  --pore-activity FILE [FILE ...]   One or more pore_activity_*.csv files.
                                    If omitted, auto-discovered in CWD.
  --throughput   FILE [FILE ...]    One or more throughput_*.csv files.
                                    If omitted, auto-discovered in CWD.
  --labels       LABEL [LABEL ...]  Human-readable sample labels (must match
                                    the file order). Auto-derived if omitted.
  --title        TEXT               Report title (default: "ONT Sequencing Report").
  --out          PATH               Output HTML file (default: report.html).
  --sample-hz    INT                Sampling interval for the Active-Pores plot,
                                    in minutes (default: 5).

Examples
--------
  # Auto-discover all CSVs in the current directory
  python generate_report.py

  # Explicit files with custom labels
  python generate_report.py \\
      --pore-activity pore_activity_PBC87221_*.csv pore_activity_PBE26402_*.csv \\
      --throughput    throughput_PBC87221_*.csv    throughput_PBE26402_*.csv \\
      --labels "prom1" "prom2" \\
      --out my_report.html
"""

import argparse
import glob
import json
import csv
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# COLOUR PALETTE
# ---------------------------------------------------------------------------
SAMPLE_COLORS = [
    "#E74C3C", "#3498DB", "#2ECC71", "#F39C12",
    "#9B59B6", "#1ABC9C", "#E67E22", "#34495E",
    "#E91E63", "#00BCD4",
]

STATE_COLORS = {
    "strand":                       "#4CAF50",
    "adapter":                      "#81C784",
    "pore":                         "#2196F3",
    "unblocking":                   "#FFC107",
    "unclassified":                 "#9E9E9E",
    "multiple":                     "#E91E63",
    "saturated":                    "#9C27B0",
    "no_pore":                      "#F44336",
    "unavailable":                  "#607D8B",
    "disabled":                     "#795548",
    "locked":                       "#FF5722",
    "membrane":                     "#00BCD4",
    "pending_manual_reset":         "#FFEB3B",
    "pending_mux_change":           "#FF9800",
    "unclassified_following_reset": "#BDBDBD",
    "unknown_negative":             "#F48FB1",
    "unknown_positive":             "#80DEEA",
    "zero":                         "#ECEFF1",
}

# States considered "active" (sequencing-productive or pore-available)
ACTIVE_STATES = {"strand", "adapter", "pore", "unblocking", "unclassified",
                 "multiple", "saturated"}

# States considered "blocking" (pore occupied by multiple strands or saturated)
BLOCKING_STATES = {"multiple", "saturated"}

# States where the pore is simply not available / not loaded — excluded from
# the blocking-ratio denominator so the metric reflects "of active pore time,
# what fraction was spent blocked?" rather than being diluted by idle channels.
IDLE_STATES = {
    "no_pore", "disabled", "zero", "unavailable",
    "locked", "membrane", "pending_manual_reset", "pending_mux_change",
}


# ---------------------------------------------------------------------------
# CSV PARSING
# ---------------------------------------------------------------------------

def _read_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)
    return rows


def parse_pore_activity(path):
    rows = _read_csv(path)
    time_map = {}
    for row in rows:
        state = row["Channel State"].strip()
        minute = int(row["Experiment Time (minutes)"])
        samples = int(row["State Time (samples)"])
        if minute not in time_map:
            time_map[minute] = {}
        time_map[minute][state] = time_map[minute].get(state, 0) + samples
    all_states = sorted({row["Channel State"].strip() for row in rows})
    return {"time_to_state_time": time_map, "all_states": all_states}


def parse_throughput(path):
    rows = _read_csv(path)
    if not rows:
        return {}
    last = rows[-1]
    total_reads = int(last.get("Reads", 0))
    passed_reads = int(last.get("Basecalled Reads Passed", 0))
    failed_reads = int(last.get("Basecalled Reads Failed", 0))
    total_bases = int(last.get("Basecalled Bases", 0))
    yield_gb = round(total_bases / 1e9, 2)
    return {
        "rows": rows,
        "last_row": last,
        "total_reads": total_reads,
        "passed_reads": passed_reads,
        "failed_reads": failed_reads,
        "total_bases": total_bases,
        "yield_gb": yield_gb,
    }


# ---------------------------------------------------------------------------
# METRICS COMPUTATION
# ---------------------------------------------------------------------------

def extract_flow_cell_id(pa_path):
    """Extract the flow cell ID from a pore_activity filename.

    e.g. pore_activity_PBC87221_run1.csv  ->  PBC87221
    """
    stem = Path(pa_path).stem.replace("pore_activity_", "")
    return stem.split("_")[0]


def compute_sample_metrics(label, pa_data, tp_data, color, flow_cell_id=""):
    time_map = pa_data["time_to_state_time"]
    minutes_sorted = sorted(time_map.keys())

    # Active pores per minute (proxy: active_samples / 300000 channels per min)
    active_per_min = {}
    for m in minutes_sorted:
        state_map = time_map[m]
        active_samples = sum(state_map.get(s, 0) for s in ACTIVE_STATES)
        active_per_min[m] = round(active_samples / 300_000)

    avg_active = round(sum(active_per_min.values()) / max(len(active_per_min), 1))

    def _active_at(minute):
        candidates = [m for m in minutes_sorted if m <= minute]
        return active_per_min[candidates[-1]] if candidates else 0

    pores_24h = _active_at(1440)
    pores_48h = _active_at(2880)

    # Blocking ratio — fraction of *active* pore time (i.e. excluding idle/dead
    # channels) spent in a truly blocked state (multiple strands or saturated).
    blocking_fracs = []
    for m in minutes_sorted:
        state_map = time_map[m]
        blocking = sum(state_map.get(s, 0) for s in BLOCKING_STATES)
        idle     = sum(state_map.get(s, 0) for s in IDLE_STATES)
        active_time = sum(state_map.values()) - idle
        if active_time > 0:
            blocking_fracs.append(blocking / active_time)
    blocking_pct = round(100 * sum(blocking_fracs) / max(len(blocking_fracs), 1), 2)

    # Throughput metrics
    yield_gb = tp_data.get("yield_gb", 0)
    total_reads = tp_data.get("total_reads", 0)
    passed_reads = tp_data.get("passed_reads", 0)
    pass_rate = round(100 * passed_reads / max(total_reads, 1), 2)
    yield_per_pore_mb = round(yield_gb * 1000 / max(avg_active, 1), 2)

    # Per-state average percentage across all timepoints
    state_avg_pct = {}
    for state in pa_data["all_states"]:
        fracs = []
        for m in minutes_sorted:
            state_map = time_map[m]
            total = sum(state_map.values()) or 1
            fracs.append(state_map.get(state, 0) / total * 100)
        state_avg_pct[state] = round(sum(fracs) / max(len(fracs), 1), 2)

    return {
        "label": label,
        "flow_cell_id": flow_cell_id,
        "color": color,
        "yield_gb": yield_gb,
        "total_reads": total_reads,
        "passed_reads": passed_reads,
        "pass_rate": pass_rate,
        "avg_active": avg_active,
        "yield_per_pore_mb": yield_per_pore_mb,
        "blocking_pct": blocking_pct,
        "pores_24h": pores_24h,
        "pores_48h": pores_48h,
        "active_per_min": active_per_min,
        "state_avg_pct": state_avg_pct,
        "all_states": pa_data["all_states"],
    }


# ---------------------------------------------------------------------------
# CHART DATA BUILDERS
# ---------------------------------------------------------------------------

def build_active_pores_plot_data(samples, sample_hz):
    series = []
    for s in samples:
        apm = s["active_per_min"]
        minutes = sorted(apm.keys())
        points = [{"x": m, "y": apm[m]} for m in minutes if m % sample_hz == 0]
        series.append({"name": s["label"], "color": s["color"], "data": points})
    return {"title": "Active Pores Over Time (Pore Stability)", "samples": series}


def build_cumulative_yield_data(samples, tp_data_list):
    series = []
    for s, tp in zip(samples, tp_data_list):
        rows = tp.get("rows", [])
        points = [
            {"x": int(r.get("Experiment Time (minutes)", 0)),
             "y": round(int(r.get("Basecalled Bases", 0)) / 1e9, 4)}
            for r in rows
        ]
        series.append({"name": s["label"], "color": s["color"], "data": points})
    return {"title": "Cumulative Yield Over Time (Gb)", "samples": series}


def build_avg_pore_activity_data(samples):
    seen = set()
    all_states_union = []
    for s in samples:
        for state in s["all_states"]:
            if state not in seen:
                all_states_union.append(state)
                seen.add(state)

    priority = [
        "strand", "adapter", "pore", "unblocking", "unclassified",
        "multiple", "saturated", "no_pore", "unavailable", "disabled",
        "locked", "membrane", "pending_manual_reset", "pending_mux_change",
        "unclassified_following_reset", "unknown_negative", "unknown_positive", "zero",
    ]
    ordered = [s for s in priority if s in seen] + \
              [s for s in all_states_union if s not in priority]

    sample_rows = []
    for s in samples:
        row = {"Sample": s["label"]}
        row.update({state: s["state_avg_pct"].get(state, 0.0) for state in ordered})
        sample_rows.append(row)

    colors = {state: STATE_COLORS.get(state, "#BDBDBD") for state in ordered}
    return {
        "title": "Average Pore Activity per Sample",
        "samples": sample_rows,
        "states": ordered,
        "colors": colors,
    }


# ---------------------------------------------------------------------------
# HTML TEMPLATE
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<meta name="description" content="Interactive ONT sequencing run comparison report.">
<style>
  *,*::before,*::after{{box-sizing:border-box;}}
  body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
        max-width:1060px;margin:0 auto;padding:24px 28px;color:#2d3748;
        line-height:1.55;background:#f7f8fa;}}
  h2{{color:#1a202c;border-bottom:2px solid #e2e8f0;padding-bottom:8px;margin-top:40px;}}
  h3{{color:#2d3748;margin-top:28px;font-size:1.05rem;}}
  table{{border-collapse:collapse;width:100%;margin:14px 0;font-size:13.5px;
         background:#fff;box-shadow:0 1px 4px rgba(0,0,0,.07);border-radius:6px;overflow:hidden;}}
  th{{background:#edf2f7;color:#2d3748;font-weight:600;text-align:left;
      padding:10px 12px;border-bottom:2px solid #cbd5e0;}}
  td{{padding:9px 12px;border-bottom:1px solid #e9ecef;}}
  tr:last-child td{{border-bottom:none;}}
  tr:hover td{{background:#f0f4ff;}}
  .plot-container{{margin:20px 0;background:#fff;border-radius:8px;
                   box-shadow:0 1px 4px rgba(0,0,0,.08);border:1px solid #e2e8f0;overflow:hidden;}}
  .footer{{margin-top:48px;border-top:1px solid #e2e8f0;padding-top:14px;
           color:#a0aec0;font-size:12px;text-align:center;}}
</style>
</head>
<body>

<h2>{title}</h2>

<h3>Metrics Summary &amp; Pore Retention</h3>
<table>
  <thead>
    <tr>
      <th>Sample</th>
      <th>Flow Cell ID</th>
      <th style="text-align:right">Yield (Gb)</th>
      <th style="text-align:right">Total Reads</th>
      <th style="text-align:right">Passed Reads</th>
      <th style="text-align:right">Pass Rate (%)</th>
      <th style="text-align:right">Avg Active Pores</th>
      <th style="text-align:right">Yield / Pore (Mb)</th>
      <th style="text-align:right">Blocking (%)</th>
      <th style="text-align:right">Pores 24 h</th>
      <th style="text-align:right">Pores 48 h</th>
    </tr>
  </thead>
  <tbody>
{summary_rows}
  </tbody>
</table>

<h3>Active Pores Over Time (Pore Stability)</h3>
<div class="plot-container">
  <iframe srcdoc="{active_pores_iframe}"
          width="100%" height="680" frameborder="0"
          scrolling="no" style="border:none;display:block;"></iframe>
</div>

<h3>Cumulative Yield Over Time</h3>
<div class="plot-container">
  <iframe srcdoc="{yield_iframe}"
          width="100%" height="480" frameborder="0"
          scrolling="no" style="border:none;display:block;"></iframe>
</div>

<h3>Average Pore Activity</h3>
<div class="plot-container">
  <iframe srcdoc="{pore_activity_iframe}"
          width="100%" height="820" frameborder="0"
          scrolling="no" style="border:none;display:block;"></iframe>
</div>

<div class="footer">
  Generated on {generated_on} &mdash; {title} &mdash;
  <a href="https://github.com/angelovangel/etc" target="_blank" rel="noopener"
     style="color:#64748b;text-decoration:none;">github.com/angelovangel/etc</a>
</div>
</body>
</html>
"""


LINE_CHART_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{chart_title}</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:#fff;margin:0;padding:20px;color:#333;}}
    #chart-container{{max-width:960px;margin:0 auto;position:relative;}}
    .line{{fill:none;stroke-width:2.5px;opacity:0.85;}}
    .tooltip{{position:absolute;padding:10px;font-size:13px;
              background:rgba(255,255,255,.95);border:1px solid #ddd;border-radius:4px;
              pointer-events:none;opacity:0;box-shadow:0 4px 6px rgba(0,0,0,.1);
              white-space:nowrap;transition:opacity .2s;z-index:10;}}
    .grid line{{stroke:#ddd;stroke-dasharray:4,4;}}
    .grid path{{stroke-width:0;}}
    .axis text{{font-size:12px;fill:#555;}}
    .axis path,.axis line{{stroke:#bbb;}}
    .hover-line{{stroke:#888;stroke-dasharray:4,4;stroke-width:1px;pointer-events:none;opacity:0;}}
    h2{{text-align:center;color:#222;margin-bottom:6px;}}
    .legend-item{{cursor:pointer;}}
    .legend-text{{font-size:12px;fill:#333;user-select:none;}}
    .toggle-bar{{display:flex;justify-content:center;gap:0;margin-bottom:10px;}}
    .toggle-btn{{
      padding:5px 18px;font-size:12px;cursor:pointer;border:1px solid #bbb;
      background:#f4f4f4;color:#444;transition:all .2s;user-select:none;
    }}
    .toggle-btn:first-child{{border-radius:6px 0 0 6px;}}
    .toggle-btn:last-child{{border-radius:0 6px 6px 0;}}
    .toggle-btn.active{{background:#2563eb;color:#fff;border-color:#2563eb;font-weight:600;}}
  </style>
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title">{chart_title}</h2>
    <div class="toggle-bar">
      <button class="toggle-btn active" id="btn-yield" onclick="setMode('yield')">Yield Over Time</button>
      <button class="toggle-btn" id="btn-rate" onclick="setMode('rate')">Rate (Gb/min)</button>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = {plot_data_json};

    // First derivative of cumulative yield: windowed rate in Gb/min.
    // Linearly interpolates the cumulative curve at regular BIN_MIN intervals,
    // then computes slope = ΔY / BIN_MIN — giving a smooth, readable signal
    // instead of noisy point-by-point differences on dense throughput data.
    const BIN_MIN = 100;
    function interpY(data, t) {{
      if (t <= data[0].x) return data[0].y;
      if (t >= data[data.length-1].x) return data[data.length-1].y;
      let lo = 0, hi = data.length - 1;
      while (hi - lo > 1) {{
        const mid = (lo + hi) >> 1;
        if (data[mid].x <= t) lo = mid; else hi = mid;
      }}
      const frac = (t - data[lo].x) / (data[hi].x - data[lo].x);
      return data[lo].y + frac * (data[hi].y - data[lo].y);
    }}
    plotData.samples.forEach(s => {{
      if (s.data.length < 2) {{ s.rateData = []; return; }}
      const xMax = s.data[s.data.length - 1].x;
      const pts = [];
      for (let t = BIN_MIN; t <= xMax + BIN_MIN / 2; t += BIN_MIN) {{
        const t1 = Math.min(t, xMax);
        const t0 = t1 - BIN_MIN;
        const rate = (interpY(s.data, t1) - interpY(s.data, t0)) / BIN_MIN;
        pts.push({{x: t1, y: Math.max(0, rate)}});
      }}
      s.rateData = pts;
    }});

    let mode = 'yield'; // 'yield' | 'rate'

    const margin = {{top:20,right:180,bottom:50,left:65}};
    const W = 960 - margin.left - margin.right;
    const H = {chart_height} - margin.top - margin.bottom;

    const svg = d3.select("#chart").append("svg")
        .attr("width",  W + margin.left + margin.right)
        .attr("height", H + margin.top  + margin.bottom)
      .append("g")
        .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const tooltip = d3.select("#tooltip");

    const allPoints = plotData.samples.flatMap(s => s.data);
    const xDomain = d3.extent(allPoints, d => d.x);
    const x = d3.scaleLinear().domain(xDomain).range([0, W]);

    // Axes groups — reused and updated on each redraw
    const gGrid  = svg.append("g").attr("class","grid");
    const gXAxis = svg.append("g").attr("class","axis")
        .attr("transform", `translate(0,${{H}})`);
    const gYAxis = svg.append("g").attr("class","axis");

    svg.append("text")
       .attr("x", W/2).attr("y", H+42)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555")
       .text("Experiment Time (minutes)");

    const yLabel = svg.append("text")
       .attr("transform","rotate(-90)")
       .attr("x", -H/2).attr("y", -55)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555");

    const lineGen = d3.line().x(d=>x(d.x)).curve(d3.curveMonotoneX);

    // Create one path + dots group per sample (updated on redraw)
    plotData.samples.forEach(sample => {{
      const g = svg.append("g")
          .attr("class", `sample-group sample-${{sample.name.replace(/\W/g,"_")}}`);
      g.append("path").attr("class","line").attr("stroke", sample.color);
      g.append("g").attr("class","dots");
    }});

    const hoverLine = svg.append("line")
        .attr("class","hover-line").attr("y1",0).attr("y2",H);

    const overlay = svg.append("rect")
        .attr("width",W).attr("height",H).attr("fill","none").attr("pointer-events","all");

    // Legend (static — click to show/hide, hover to highlight)
    const legend = svg.append("g").attr("transform",`translate(${{W+14}},0)`);
    plotData.samples.forEach((s,i) => {{
      const row = legend.append("g").attr("class","legend-item")
          .attr("transform",`translate(0,${{i*24}})`)
          .on("mouseover", () => {{
            svg.selectAll(".sample-group").style("opacity",0.1);
            svg.select(`.sample-${{s.name.replace(/\W/g,"_")}}`).style("opacity",1).raise();
          }})
          .on("mouseout", () => svg.selectAll(".sample-group").style("opacity",1))
          .on("click", function() {{
            s.hidden = !s.hidden;
            svg.select(`.sample-${{s.name.replace(/\W/g,"_")}}`)
               .style("display", s.hidden ? "none" : null);
            d3.select(this).select("circle").style("fill", s.hidden ? "#ccc" : s.color);
            d3.select(this).select("text")
              .style("text-decoration", s.hidden ? "line-through" : "none")
              .style("fill", s.hidden ? "#999" : "#333");
          }});
      row.append("circle").attr("r",5).style("fill",s.color);
      row.append("text").attr("class","legend-text").attr("x",10).attr("y",4).text(s.name);
    }});

    let y = d3.scaleLinear().range([H, 0]);

    function draw() {{
      const isRate = mode === 'rate';
      const activeKey = isRate ? 'rateData' : 'data';
      const label = isRate ? 'Rate (Gb/min)' : '{y_label}';
      document.getElementById("chart-title").textContent =
        isRate ? 'Yield Rate Over Time (Gb/min)' : '{chart_title}';

      const activePoints = plotData.samples.flatMap(s => s[activeKey]);
      const yMax = d3.max(activePoints, d => d.y);
      const yMin = isRate ? d3.min(activePoints, d => d.y) : 0;
      y.domain([Math.min(0, yMin), yMax * 1.05]);

      gGrid.call(d3.axisLeft(y).tickSize(-W).tickFormat(""));
      gXAxis.call(d3.axisBottom(x).ticks(10).tickFormat(d => d + " min"));
      gYAxis.call(d3.axisLeft(y));
      yLabel.text(label);

      lineGen.y(d => y(d.y));

      plotData.samples.forEach(sample => {{
        const cls = `.sample-${{sample.name.replace(/\W/g,"_")}}`;
        svg.select(cls + " path")
           .datum(sample[activeKey])
           .attr("d", lineGen);
        // Dot markers — shown only in rate mode
        svg.select(cls + " .dots").selectAll("circle")
           .data(isRate ? sample[activeKey] : [])
           .join("circle")
           .attr("r", 3.5)
           .attr("cx", d => x(d.x))
           .attr("cy", d => y(d.y))
           .attr("fill", sample.color)
           .attr("stroke", "#fff")
           .attr("stroke-width", 1.5);
      }});

      overlay
        .on("mousemove", function(event) {{
          const [mx] = d3.pointer(event);
          const xVal = x.invert(mx);
          hoverLine.attr("x1",mx).attr("x2",mx).style("opacity",1);
          let html = `<strong>t = ${{Math.round(xVal)}} min</strong><br/>`;
          plotData.samples.forEach(s => {{
            const bisect = d3.bisector(d=>d.x).left;
            const data = s[activeKey];
            const idx = bisect(data, xVal);
            const pt = data[idx] || data[data.length-1];
            if (pt) html += `<span style="color:${{s.color}}">&#9632;</span> ${{s.name}}: <b>${{pt.y.toLocaleString()}}</b><br/>`;
          }});
          tooltip.transition().duration(40).style("opacity",1);
          tooltip.html(html)
                 .style("left",(event.pageX+15)+"px")
                 .style("top",Math.min(event.pageY-28, window.innerHeight-160)+"px");
        }})
        .on("mouseout", () => {{
          hoverLine.style("opacity",0);
          tooltip.transition().duration(200).style("opacity",0);
        }});
    }}

    function setMode(m) {{
      mode = m;
      document.getElementById("btn-yield").classList.toggle("active", m === 'yield');
      document.getElementById("btn-rate").classList.toggle("active",  m === 'rate');
      draw();
    }}

    draw();
  </script>
</body>
</html>"""


ACTIVE_PORES_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{chart_title}</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:#fff;margin:0;padding:20px;color:#333;}}
    #chart-container{{max-width:960px;margin:0 auto;position:relative;}}
    .line{{fill:none;stroke-width:2.5px;opacity:0.85;}}
    .tooltip{{position:absolute;padding:10px;font-size:13px;
              background:rgba(255,255,255,.95);border:1px solid #ddd;border-radius:4px;
              pointer-events:none;opacity:0;box-shadow:0 4px 6px rgba(0,0,0,.1);
              white-space:nowrap;transition:opacity .2s;z-index:10;}}
    .grid line{{stroke:#ddd;stroke-dasharray:4,4;}}
    .grid path{{stroke-width:0;}}
    .axis text{{font-size:12px;fill:#555;}}
    .axis path,.axis line{{stroke:#bbb;}}
    .hover-line{{stroke:#888;stroke-dasharray:4,4;stroke-width:1px;pointer-events:none;opacity:0;}}
    h2{{text-align:center;color:#222;margin-bottom:6px;}}
    .legend-item{{cursor:pointer;}}
    .legend-text{{font-size:12px;fill:#333;user-select:none;}}
    .toggle-bar{{display:flex;justify-content:center;gap:0;margin-bottom:10px;}}
    .toggle-btn{{
      padding:5px 18px;font-size:12px;cursor:pointer;border:1px solid #bbb;
      background:#f4f4f4;color:#444;transition:all .2s;user-select:none;
    }}
    .toggle-btn:first-child{{border-radius:6px 0 0 6px;}}
    .toggle-btn:last-child{{border-radius:0 6px 6px 0;}}
    .toggle-btn.active{{background:#2563eb;color:#fff;border-color:#2563eb;font-weight:600;}}
  </style>
</head>
<body>
  <div id="chart-container">
    <h2>{chart_title}</h2>
    <div class="toggle-bar">
      <button class="toggle-btn active" id="btn-abs" onclick="setMode('abs')">Absolute</button>
      <button class="toggle-btn" id="btn-rel" onclick="setMode('rel')">Relative (%)</button>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = {plot_data_json};

    // Pre-compute baseline (max y value) for each sample — 100% = peak pore count
    plotData.samples.forEach(s => {{
      s.baseline = d3.max(s.data, d => d.y) || 1;
      s.relData = s.data.map(d => ({{x: d.x, y: +(d.y / s.baseline * 100).toFixed(2)}}));
    }});

    let mode = 'abs'; // 'abs' | 'rel'

    const margin = {{top:20,right:180,bottom:50,left:65}};
    const W = 960 - margin.left - margin.right;
    const H = {chart_height} - margin.top - margin.bottom;

    const svg = d3.select("#chart").append("svg")
        .attr("width",  W + margin.left + margin.right)
        .attr("height", H + margin.top  + margin.bottom)
      .append("g")
        .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const tooltip = d3.select("#tooltip");
    const allPoints = plotData.samples.flatMap(s => s.data);
    const xDomain = d3.extent(allPoints, d => d.x);
    const x = d3.scaleLinear().domain(xDomain).range([0, W]);
    const y = d3.scaleLinear().range([H, 0]);

    const gridG  = svg.append("g").attr("class","grid");
    const xAxisG = svg.append("g").attr("class","axis").attr("transform",`translate(0,${{H}})`);
    const yAxisG = svg.append("g").attr("class","axis");

    xAxisG.call(d3.axisBottom(x).ticks(10).tickFormat(d => d + " min"));

    svg.append("text")
       .attr("x", W/2).attr("y", H+42)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555")
       .text("Experiment Time (minutes)");

    const yLabel = svg.append("text")
       .attr("transform","rotate(-90)")
       .attr("x", -H/2).attr("y", -55)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555");

    const lineGen = d3.line().x(d=>x(d.x)).y(d=>y(d.y)).curve(d3.curveMonotoneX);

    plotData.samples.forEach(sample => {{
      const g = svg.append("g")
          .attr("class", `sample-group sample-${{sample.name.replace(/\W/g,"_")}}`);
      g.append("path").attr("class","line").attr("stroke", sample.color);
    }});

    const hoverLine = svg.append("line")
        .attr("class","hover-line").attr("y1",0).attr("y2",H);

    svg.append("rect")
        .attr("width",W).attr("height",H).attr("fill","none").attr("pointer-events","all")
        .on("mousemove", function(event) {{
          const [mx] = d3.pointer(event);
          const xVal = x.invert(mx);
          hoverLine.attr("x1",mx).attr("x2",mx).style("opacity",1);
          const isRel = mode === 'rel';
          let html = `<strong>t = ${{Math.round(xVal)}} min</strong><br/>`;
          plotData.samples.forEach(s => {{
            const pts = isRel ? s.relData : s.data;
            const bisect = d3.bisector(d=>d.x).left;
            const idx = bisect(pts, xVal);
            const d = pts[idx] || pts[pts.length-1];
            if (d) {{
              const val = isRel ? d.y.toFixed(1) + '%' : d.y.toLocaleString();
              html += `<span style="color:${{s.color}}">&#9632;</span> ${{s.name}}: <b>${{val}}</b><br/>`;
            }}
          }});
          tooltip.transition().duration(40).style("opacity",1);
          tooltip.html(html)
                 .style("left",(event.pageX+15)+"px")
                 .style("top",Math.min(event.pageY-28, window.innerHeight-160)+"px");
        }})
        .on("mouseout", () => {{
          hoverLine.style("opacity",0);
          tooltip.transition().duration(200).style("opacity",0);
        }});

    const legend = svg.append("g").attr("transform",`translate(${{W+14}},0)`);
    plotData.samples.forEach((s,i) => {{
      const row = legend.append("g").attr("class","legend-item")
          .attr("transform",`translate(0,${{i*24}})`)
          .on("mouseover", () => {{
            svg.selectAll(".sample-group").style("opacity",0.1);
            svg.select(`.sample-${{s.name.replace(/\W/g,"_")}}`).style("opacity",1).raise();
          }})
          .on("mouseout", () => svg.selectAll(".sample-group").style("opacity",1))
          .on("click", function() {{
            s.hidden = !s.hidden;
            svg.select(`.sample-${{s.name.replace(/\W/g,"_")}}`)
               .style("display", s.hidden ? "none" : null);
            d3.select(this).select("circle").style("fill", s.hidden ? "#ccc" : s.color);
            d3.select(this).select("text")
              .style("text-decoration", s.hidden ? "line-through" : "none")
              .style("fill", s.hidden ? "#999" : "#333");
          }});
      row.append("circle").attr("r",5).style("fill",s.color);
      row.append("text").attr("class","legend-text").attr("x",10).attr("y",4).text(s.name);
    }});

    function setMode(m) {{
      mode = m;
      document.getElementById('btn-abs').classList.toggle('active', m==='abs');
      document.getElementById('btn-rel').classList.toggle('active', m==='rel');
      render();
    }}

    function render() {{
      const isRel = mode === 'rel';
      const pts = plotData.samples.flatMap(s => isRel ? s.relData : s.data);
      const yMax = d3.max(pts, d => d.y);
      y.domain([0, yMax * 1.05]);

      const fmt = isRel ? d => d + '%' : d => d.toLocaleString();
      yAxisG.transition().duration(400).call(d3.axisLeft(y).tickFormat(fmt));
      gridG.transition().duration(400)
           .call(d3.axisLeft(y).tickSize(-W).tickFormat(""));
      yLabel.text(isRel ? 'Active Pores (% of start)' : 'Active Pores');

      plotData.samples.forEach(sample => {{
        const data = isRel ? sample.relData : sample.data;
        svg.select(`.sample-${{sample.name.replace(/\W/g,"_")}} path`)
           .datum(data)
           .transition().duration(400)
           .attr("d", lineGen);
      }});
    }}

    render();
  </script>
</body>
</html>"""


STACKED_BAR_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Average Pore Activity</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:#fff;margin:0;padding:20px;color:#333;}}
    #chart-container{{max-width:960px;margin:0 auto;position:relative;}}
    .tooltip{{position:absolute;padding:10px;font-size:13px;background:rgba(255,255,255,.95);
              border:1px solid #ddd;border-radius:4px;pointer-events:none;opacity:0;
              box-shadow:0 4px 6px rgba(0,0,0,.1);white-space:nowrap;transition:opacity .2s;z-index:10;}}
    .axis text{{font-size:12px;fill:#555;}}
    .axis path,.axis line{{stroke:#bbb;}}
    .grid line{{stroke:#eee;}}
    h2{{text-align:center;color:#222;margin-bottom:4px;}}
    p.subtitle{{text-align:center;color:#666;font-size:13px;margin-bottom:20px;}}
    .legend-container{{display:flex;flex-wrap:wrap;justify-content:center;gap:8px;
                       margin-top:18px;padding:12px;background:#fcfcfc;border-radius:8px;
                       border:1px solid #eee;}}
    .legend-item{{display:flex;align-items:center;font-size:12px;color:#444;cursor:pointer;
                  padding:5px 10px;border-radius:5px;transition:all .2s;user-select:none;
                  border:1px solid transparent;}}
    .legend-item:hover{{background:#f0f0f0;border-color:#ddd;}}
    .legend-item.hidden{{opacity:.35;}}
    .legend-item.hidden .legend-text{{text-decoration:line-through;color:#888;}}
    .legend-color{{width:13px;height:13px;margin-right:7px;border-radius:3px;
                   border:1px solid rgba(0,0,0,.12);}}
  </style>
</head>
<body>
  <div id="chart-container">
    <h2>Average Pore Activity per Sample</h2>
    <p class="subtitle">Click legend items to toggle states</p>
    <div id="chart"></div>
    <div class="legend-container" id="legend"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = {plot_data_json};

    const margin = {{top:20,right:30,bottom:100,left:65}};
    const W = 960 - margin.left - margin.right;
    const H = 500 - margin.top - margin.bottom;

    const svg = d3.select("#chart").append("svg")
        .attr("width",  W + margin.left + margin.right)
        .attr("height", H + margin.top  + margin.bottom)
      .append("g")
        .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const tooltip = d3.select("#tooltip");

    const x = d3.scaleBand().domain(plotData.samples.map(d=>d.Sample)).range([0,W]).padding(0.3);
    const y = d3.scaleLinear().range([H,0]);
    const stack = d3.stack();
    let visibleStates = [...plotData.states];

    const xAxis = svg.append("g").attr("class","axis").attr("transform",`translate(0,${{H}})`);
    xAxis.call(d3.axisBottom(x)).selectAll("text")
         .attr("transform","rotate(-30)").style("text-anchor","end").attr("dy","0.8em");

    const gridGroup = svg.append("g").attr("class","grid");
    const yAxisGroup = svg.append("g").attr("class","axis");

    svg.append("text").attr("transform","rotate(-90)")
       .attr("x",-H/2).attr("y",-55)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555")
       .text("% Time in State");

    function updateChart() {{
      const totals = plotData.samples.map(d => visibleStates.reduce((s,k)=>s+(d[k]||0),0));
      const maxTotal = d3.max(totals) || 100;
      y.domain([0, Math.max(100, maxTotal*1.05)]);
      yAxisGroup.transition().duration(400).call(d3.axisLeft(y).tickFormat(d=>d+"%"));
      gridGroup.transition().duration(400).call(d3.axisLeft(y).tickSize(-W).tickFormat(""));

      stack.keys(visibleStates);
      const stackedData = stack(plotData.samples);
      const layers = svg.selectAll(".layer").data(stackedData, d=>d.key);
      layers.exit().remove();
      const layersE = layers.enter().append("g").attr("class","layer")
            .attr("fill", d=>plotData.colors[d.key]||"#ccc");
      const allLayers = layersE.merge(layers);

      const rects = allLayers.selectAll("rect").data(d=>d, d=>d.data.Sample);
      rects.exit().remove();
      rects.enter().append("rect")
           .attr("x", d=>x(d.data.Sample)).attr("width",x.bandwidth())
           .attr("y",H).attr("height",0)
           .merge(rects)
           .on("mouseover", function(event,d) {{
             const key = d3.select(this.parentNode).datum().key;
             const color = plotData.colors[key] || "#ccc";
             let html = `<strong>${{key}}</strong><br/>`;
             plotData.samples.forEach(s => {{
               const val = (s[key] || 0).toFixed(2);
               html += `<span style="display:inline-block;width:10px;height:10px;`
                     + `background:${{color}};border-radius:2px;margin-right:6px;`
                     + `vertical-align:middle;border:1px solid rgba(0,0,0,.15);"></span>`
                     + `${{s.Sample}}: <b>${{val}}%</b><br/>`;
             }});
             tooltip.transition().duration(40).style("opacity",1);
             tooltip.html(html)
                    .style("left",(event.pageX+10)+"px")
                    .style("top",(event.pageY-10)+"px");
             d3.select(this).attr("opacity",0.7);
           }})
           .on("mousemove", function(event) {{
             tooltip.style("left",(event.pageX+10)+"px")
                    .style("top",(event.pageY-10)+"px");
           }})
           .on("mouseout", function() {{
             tooltip.transition().duration(200).style("opacity",0);
             d3.select(this).attr("opacity",1);
           }})
           .transition().duration(400)
           .attr("y", d=>y(d[1])).attr("height", d=>Math.max(0,y(d[0])-y(d[1])));
    }}

    updateChart();

    const legendContainer = d3.select("#legend");
    plotData.states.forEach(state => {{
      const item = legendContainer.append("div").attr("class","legend-item")
          .on("click", function() {{
            const idx = visibleStates.indexOf(state);
            if (idx > -1) {{
              if (visibleStates.length > 1) {{
                visibleStates.splice(idx,1);
                d3.select(this).classed("hidden",true);
              }}
            }} else {{
              visibleStates.push(state);
              visibleStates.sort((a,b)=>plotData.states.indexOf(a)-plotData.states.indexOf(b));
              d3.select(this).classed("hidden",false);
            }}
            updateChart();
          }});
      item.append("div").attr("class","legend-color").style("background-color",plotData.colors[state]||"#ccc");
      item.append("span").attr("class","legend-text").text(state);
    }});
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# RENDERING HELPERS
# ---------------------------------------------------------------------------

def _html_escape_attr(s):
    return (s.replace("&", "&amp;")
             .replace('"', "&quot;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace("'", "&#x27;"))


def _render_line_chart_iframe(plot_data, chart_height, y_label):
    title = plot_data.get("title", "")
    json_str = json.dumps(plot_data, separators=(",", ":"))
    html = LINE_CHART_IFRAME.format(
        chart_title=title,
        plot_data_json=json_str,
        chart_height=chart_height,
        y_label=y_label,
    )
    return _html_escape_attr(html)


def _render_active_pores_iframe(plot_data, chart_height):
    title = plot_data.get("title", "")
    json_str = json.dumps(plot_data, separators=(",", ":"))
    html = ACTIVE_PORES_IFRAME.format(
        chart_title=title,
        plot_data_json=json_str,
        chart_height=chart_height,
    )
    return _html_escape_attr(html)


def _render_stacked_bar_iframe(plot_data):
    json_str = json.dumps(plot_data, separators=(",", ":"))
    html = STACKED_BAR_IFRAME.format(plot_data_json=json_str)
    return _html_escape_attr(html)


def _render_summary_row(s):
    return (
        f"    <tr>"
        f"<td><strong>{s['label']}</strong></td>"
        f"<td>{s.get('flow_cell_id', '')}</td>"
        f"<td style='text-align:right'>{s['yield_gb']:.2f}</td>"
        f"<td style='text-align:right'>{s['total_reads']:,}</td>"
        f"<td style='text-align:right'>{s['passed_reads']:,}</td>"
        f"<td style='text-align:right'>{s['pass_rate']:.2f}%</td>"
        f"<td style='text-align:right'>{s['avg_active']:,}</td>"
        f"<td style='text-align:right'>{s['yield_per_pore_mb']:.2f}</td>"
        f"<td style='text-align:right'>{s['blocking_pct']:.2f}%</td>"
        f"<td style='text-align:right'>{s['pores_24h']:,}</td>"
        f"<td style='text-align:right'>{s['pores_48h']:,}</td>"
        f"</tr>"
    )


# ---------------------------------------------------------------------------
# FILE PAIRING
# ---------------------------------------------------------------------------

def match_pairs(pa_files, tp_files):
    def run_id(path):
        stem = Path(path).stem
        for prefix in ("pore_activity_", "throughput_"):
            stem = stem.replace(prefix, "")
        return stem

    pa_by_id = {run_id(f): f for f in pa_files}
    tp_by_id = {run_id(f): f for f in tp_files}
    shared = {run_id(f) for f in pa_files} & set(tp_by_id)
    if shared:
        # Preserve the original pa_files order so --labels align correctly.
        return [(pa_by_id[rid], tp_by_id[rid])
                for rid in (run_id(f) for f in pa_files)
                if rid in shared]
    return list(zip(pa_files, tp_files))


def derive_label(pa_path):
    stem = Path(pa_path).stem.replace("pore_activity_", "")
    return stem.split("_")[0]


# ---------------------------------------------------------------------------
# MAIN GENERATION FUNCTION
# ---------------------------------------------------------------------------

def generate_report(pa_files, tp_files, labels, title, out_path, sample_hz):
    pairs = match_pairs(pa_files, tp_files)
    if not pairs:
        sys.exit("ERROR: No matching pore_activity/throughput file pairs found.")

    print(f"Found {len(pairs)} sample pair(s).")

    samples_metrics = []
    tp_data_list = []

    for i, (pa_path, tp_path) in enumerate(pairs):
        label = labels[i] if i < len(labels) else derive_label(pa_path)
        color = SAMPLE_COLORS[i % len(SAMPLE_COLORS)]
        flow_cell_id = extract_flow_cell_id(pa_path)
        print(f"  [{i+1}] {label}")
        print(f"       pore_activity : {pa_path}")
        print(f"       throughput    : {tp_path}")

        pa_data = parse_pore_activity(pa_path)
        tp_data = parse_throughput(tp_path)
        tp_data_list.append(tp_data)

        metrics = compute_sample_metrics(label, pa_data, tp_data, color, flow_cell_id)
        samples_metrics.append(metrics)

    # Build JSON payloads
    active_pores_data = build_active_pores_plot_data(samples_metrics, sample_hz)
    yield_data        = build_cumulative_yield_data(samples_metrics, tp_data_list)
    pore_act_data     = build_avg_pore_activity_data(samples_metrics)

    # Render iframes (escaped for HTML attributes)
    active_pores_iframe  = _render_active_pores_iframe(active_pores_data, 600)
    yield_iframe         = _render_line_chart_iframe(yield_data, 400, "Yield (Gb)")
    pore_activity_iframe = _render_stacked_bar_iframe(pore_act_data)

    summary_rows = "\n".join(_render_summary_row(s) for s in samples_metrics)

    report = HTML_TEMPLATE.format(
        title=title,
        summary_rows=summary_rows,
        active_pores_iframe=active_pores_iframe,
        yield_iframe=yield_iframe,
        pore_activity_iframe=pore_activity_iframe,
        generated_on=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )

    Path(out_path).write_text(report, encoding="utf-8")
    print(f"\nReport written to: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive ONT sequencing HTML report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--pore-activity", nargs="+", metavar="FILE",
                        help="pore_activity_*.csv files (auto-discovered if omitted)")
    parser.add_argument("--throughput", nargs="+", metavar="FILE",
                        help="throughput_*.csv files (auto-discovered if omitted)")
    parser.add_argument("--labels", nargs="+", metavar="LABEL",
                        help="Human-readable sample labels (must match file order)")
    parser.add_argument("--title", default="ONT Sequencing Report",
                        help="Report title (default: 'ONT Sequencing Report')")
    parser.add_argument("--out", default="report.html",
                        help="Output HTML file (default: report.html)")
    parser.add_argument("--sample-hz", type=int, default=5, metavar="INT",
                        help="Sampling interval for Active-Pores plot in minutes (default: 5)")
    args = parser.parse_args()

    pa_files = sorted(args.pore_activity or glob.glob("pore_activity_*.csv"))
    tp_files = sorted(args.throughput   or glob.glob("throughput_*.csv"))

    if not pa_files:
        sys.exit("ERROR: No pore_activity_*.csv files found. "
                 "Use --pore-activity to specify them explicitly.")
    if not tp_files:
        sys.exit("ERROR: No throughput_*.csv files found. "
                 "Use --throughput to specify them explicitly.")

    generate_report(
        pa_files=pa_files,
        tp_files=tp_files,
        labels=args.labels or [],
        title=args.title,
        out_path=args.out,
        sample_hz=args.sample_hz,
    )


if __name__ == "__main__":
    main()
