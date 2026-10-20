#!/usr/bin/env python3
"""
generate-ont-report.py
======================
Generate an interactive HTML report for ONT (PromethION/MinION) sequencing
runs from MinKNOW's pore_activity_*.csv and throughput_*.csv files, plus
optional sequencing_summary_*.txt files for a read-length Nx plot.

Usage
-----
  generate-ont-report.py [INPUT ...] [-o OUT] [-t TITLE] [-l LABEL ...]

INPUT is any mix of MinKNOW files and run folders (searched recursively);
with no INPUT the current folder is searched. Files are matched by the run
ID in their names (pore_activity_<run>.csv, throughput_<run>.csv,
sequencing_summary_<run>.txt), so their order doesn't matter. The Nx plot
is included when every run has a sequencing summary.

Options
-------
  -o, --out PATH        Output HTML file (default: ont-report.html).
  -t, --title TEXT      Report title (default: "ONT Sequencing Report").
  -l, --labels LABEL …  Sample labels. Either NAME=LABEL, where NAME is a flow
                        cell or run ID (e.g. PBC87221=Liver), or plain labels
                        in the order the runs are listed. Default: flow cell ID.
  --sample-hz MIN       Sampling interval of the active-pores plot (default: 5).

  --pore-activity, --throughput, --sequencing-summary FILE …
                        Give files of one kind explicitly, e.g. when they were
                        renamed. Files whose names share no run ID are paired
                        in the order given.
  --nx-values JSON      Precomputed Nx values, one array per run (used by the app).

Examples
--------
  # Two MinKNOW run folders, labelled by flow cell
  generate-ont-report.py /data/run1 /data/run2 -l PBC87221=Liver PBE26402=Blood

  # All files in the current folder
  generate-ont-report.py -t "Week 42" -o week42.html
"""

import argparse
import html as html_lib
import json
import csv
import math
import statistics
import sys
from datetime import datetime
from pathlib import Path

from report_common import (D3_SOURCE_TAG, PAGE_TIP, TABLE_CSV_JS, fill, footer_html,
                           load_d3, render_chart_iframe, tip_attr, tip_html)

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

# States where a pore is able to sequence (in a strand, capturing, or open and
# waiting) — the denominator for per-pore yield rate.
SEQUENCING_STATES = {"strand", "adapter", "pore"}

# Sequencing-capable states in which the pore is occupied — numerator of pore
# occupancy, as MinKNOW defines it: (strand + adapter) / (strand + adapter + pore).
OCCUPIED_STATES = {"strand", "adapter"}

# Signal samples recorded per channel per minute at 5 kHz.
SAMPLES_PER_CHANNEL_MIN = 300_000


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


READ_LENGTH_SAMPLE_MAX = 50_000
READ_LENGTH_ESTIMATE_ROWS = 100


def parse_sequencing_summary(path):
    lengths = []
    with open(path, "rb") as fh:
        header_line = fh.readline().decode("utf-8-sig")
        header = next(csv.reader([header_line], delimiter="\t"), [])
        if "sequence_length_template" not in header:
            raise ValueError(
                f"{path}: missing required sequence_length_template column"
            )
        length_index = header.index("sequence_length_template")

        body_start = fh.tell()
        file_size = Path(path).stat().st_size
        prefix_lines = []
        for _ in range(READ_LENGTH_ESTIMATE_ROWS):
            line = fh.readline()
            if not line:
                break
            prefix_lines.append(line)

        body_size = file_size - body_start
        average_line_size = (
            sum(map(len, prefix_lines)) / len(prefix_lines)
            if prefix_lines else 1
        )
        estimated_read_count = body_size / average_line_size

        def append_length(line, line_number):
            row = next(csv.reader([line.decode("utf-8")], delimiter="\t"), [])
            raw_length = row[length_index].strip() if length_index < len(row) else ""
            try:
                length = int(raw_length)
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid sequence_length_template "
                    f"value {raw_length!r}"
                ) from exc
            if length < 0:
                raise ValueError(
                    f"{path}:{line_number}: read length cannot be negative"
                )
            if length > 0:
                lengths.append(length)

        if estimated_read_count <= READ_LENGTH_SAMPLE_MAX:
            for line_number, line in enumerate(prefix_lines, start=2):
                append_length(line, line_number)
            line_number = len(prefix_lines) + 2
            for line in fh:
                append_length(line, line_number)
                line_number += 1
        else:
            sample_count = READ_LENGTH_SAMPLE_MAX
            sample_spacing = body_size / sample_count
            for sample_index in range(sample_count):
                target = body_start + int(sample_index * sample_spacing)
                fh.seek(target)
                if target > body_start:
                    fh.seek(target - 1)
                    fh.readline()
                line = fh.readline()
                if line:
                    estimated_line_number = (
                        int((target - body_start) / average_line_size) + 2
                    )
                    append_length(line, estimated_line_number)

    if not lengths:
        raise ValueError(f"{path}: no positive read lengths found")
    return lengths


# ---------------------------------------------------------------------------
# METRICS COMPUTATION
# ---------------------------------------------------------------------------

def compute_sample_metrics(label, pa_data, tp_data, color, flow_cell_id=""):
    time_map = pa_data["time_to_state_time"]
    minutes_sorted = sorted(time_map.keys())

    # Minutes recording much less channel time than usual are incomplete (e.g.
    # the final minute when the run stopped) and would understate pore counts.
    totals = sorted(sum(time_map[m].values()) for m in minutes_sorted)
    full_minute = totals[len(totals) // 2] if totals else 0

    # Active pores per minute (proxy: active_samples / 300000 channels per min)
    active_per_min = {}
    sequencing_samples_per_min = {}
    occupied_samples_per_min = {}
    for m in minutes_sorted:
        state_map = time_map[m]
        active_samples = sum(state_map.get(s, 0) for s in ACTIVE_STATES)
        active_per_min[m] = round(active_samples / SAMPLES_PER_CHANNEL_MIN)
        if sum(state_map.values()) >= 0.9 * full_minute:
            sequencing_samples_per_min[m] = sum(state_map.get(s, 0) for s in SEQUENCING_STATES)
            occupied_samples_per_min[m] = sum(state_map.get(s, 0) for s in OCCUPIED_STATES)
    sequencing_per_min = {m: round(v / SAMPLES_PER_CHANNEL_MIN)
                          for m, v in sequencing_samples_per_min.items()}

    avg_active = round(sum(active_per_min.values()) / max(len(active_per_min), 1))
    avg_sequencing = round(sum(sequencing_per_min.values()) / max(len(sequencing_per_min), 1))

    # Time-weighted over the run, so the sparse late-run minutes count for little.
    occupancy_pct = round(100 * sum(occupied_samples_per_min.values())
                          / max(sum(sequencing_samples_per_min.values()), 1), 2)

    retention = _pore_retention(active_per_min)

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
    yield_per_pore_mb = round(yield_gb * 1000 / max(avg_sequencing, 1), 2)

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
        "occupancy_pct": occupancy_pct,
        "blocking_pct": blocking_pct,
        "half_life_h": retention["half_life_h"],
        "half_life_censored": retention["censored"],
        "retention_curve": retention["curve"],
        "active_per_min": active_per_min,
        "sequencing_per_min": sequencing_per_min,
        "sequencing_samples_per_min": sequencing_samples_per_min,
        "occupied_samples_per_min": occupied_samples_per_min,
        "state_avg_pct": state_avg_pct,
        "all_states": pa_data["all_states"],
    }


# Retention smoothing: a rolling median over this many minutes removes the
# short pore-count dips of mux scans and flushes, which would otherwise end the
# half-life early.
RETENTION_WINDOW_MIN = 61
# Resolution of the retention curve drawn in the summary table.
RETENTION_STEP_MIN = 30


def _pore_retention(active_per_min):
    """Pore half-life and a %-of-peak retention curve from active pores.

    Half-life = hours from the start until the smoothed active-pore count
    first falls below 50% of its peak (after the peak). If it never does,
    half_life_h is the run length and censored is True ("≥ N h").
    """
    minutes = sorted(active_per_min)
    if not minutes:
        return {"half_life_h": 0.0, "censored": True, "curve": []}
    values = [active_per_min[m] for m in minutes]
    half = RETENTION_WINDOW_MIN // 2
    smoothed = [statistics.median(values[max(0, i - half):i + half + 1])
                for i in range(len(values))]
    peak = max(smoothed) or 1
    peak_i = smoothed.index(max(smoothed))

    half_life_i = next((i for i in range(peak_i, len(smoothed))
                        if smoothed[i] < 0.5 * peak), None)
    censored = half_life_i is None
    half_life_min = minutes[-1] if censored else minutes[half_life_i]

    curve = [[minutes[i], round(100 * smoothed[i] / peak, 1)]
             for i in range(len(minutes))
             if minutes[i] % RETENTION_STEP_MIN == 0 or i == len(minutes) - 1]
    return {"half_life_h": round(half_life_min / 60, 1), "censored": censored, "curve": curve}


# ---------------------------------------------------------------------------
# CHART DATA BUILDERS
# ---------------------------------------------------------------------------

def build_active_pores_plot_data(samples, sample_hz):
    series = []
    for s in samples:
        apm = s["active_per_min"]
        minutes = sorted(apm.keys())
        points = [{"x": m, "y": apm[m]} for m in minutes if m % sample_hz == 0]
        series.append({"name": s["label"], "color": s["color"], "data": points,
                       "occupancy": _occupancy_series(s, minutes, sample_hz)})
    return {"title": "Active Pores Over Time (Pore Stability)", "samples": series}


def _occupancy_series(s, minutes, sample_hz):
    """Pore occupancy (%) pooled over each sample_hz-minute window."""
    seq = s["sequencing_samples_per_min"]
    occ = s["occupied_samples_per_min"]
    # As for the per-pore rate, skip windows with too few pores to be meaningful.
    min_samples = 0.05 * max(seq.values(), default=0)
    points = []
    for m in minutes:
        if m % sample_hz:
            continue
        window = [k for k in range(m - sample_hz + 1, m + 1) if k in seq]
        denom = sum(seq[k] for k in window)
        if window and denom >= min_samples * len(window) and denom > 0:
            points.append({"x": m, "y": round(100 * sum(occ[k] for k in window) / denom, 2)})
    return points


def build_cumulative_yield_data(samples, tp_data_list):
    series = []
    for s, tp in zip(samples, tp_data_list):
        rows = tp.get("rows", [])
        points = [
            {"x": int(r.get("Experiment Time (minutes)", 0)),
             "y": round(int(r.get("Basecalled Bases", 0)) / 1e9, 4)}
            for r in rows
        ]
        # Estimated Bases is derived from raw signal as it is acquired, so unlike
        # Basecalled Bases it doesn't lag behind live basecalling — used for the
        # per-pore rate, which must line up in time with pore_activity.
        estimated = [
            {"x": int(r.get("Experiment Time (minutes)", 0)),
             "y": round(int(r["Estimated Bases"]) / 1e9, 4)}
            for r in rows
        ] if rows and "Estimated Bases" in rows[0] else points
        spm = s["sequencing_per_min"]
        pores = [{"x": m, "y": spm[m]} for m in sorted(spm)]
        series.append({"name": s["label"], "color": s["color"], "data": points,
                       "estimated": estimated, "pores": pores})
    return {"title": "Cumulative Yield Over Time (Gb)", "samples": series}


def calculate_nx(lengths):
    ordered_lengths = sorted(lengths, reverse=True)
    total_bases = sum(ordered_lengths)
    if not ordered_lengths or total_bases <= 0:
        raise ValueError("Cannot calculate Nx values without positive read lengths")

    nx_values = {}
    next_percentile_index = 0
    percentiles = list(range(5, 100, 5))
    cumulative_bases = 0
    for length in ordered_lengths:
        cumulative_bases += length
        while next_percentile_index < len(percentiles):
            percentile = percentiles[next_percentile_index]
            if cumulative_bases * 100 < total_bases * percentile:
                break
            nx_values[percentile] = round(length / 1000, 3)
            next_percentile_index += 1
        if next_percentile_index == len(percentiles):
            break

    return [{"x": p, "y": nx_values[p]} for p in percentiles]


def build_nx_plot_data(samples, sequencing_summary_files):
    series = []
    for sample, summary_path in zip(samples, sequencing_summary_files):
        lengths = parse_sequencing_summary(summary_path)
        series.append({
            "name": sample["label"],
            "color": sample["color"],
            "data": calculate_nx(lengths),
        })
    return {"title": "Read Length Nx", "samples": series}


def build_nx_plot_data_from_values(samples, nx_values):
    percentiles = list(range(5, 100, 5))
    if len(nx_values) != len(samples):
        raise ValueError("Provide Nx values for every sample")

    series = []
    for sample, sample_values in zip(samples, nx_values):
        if not isinstance(sample_values, list) or len(sample_values) != len(percentiles):
            raise ValueError(
                f"{sample['label']}: expected {len(percentiles)} Nx values"
            )
        points = []
        for expected_percentile, point in zip(percentiles, sample_values):
            if not isinstance(point, dict) or point.get("x") != expected_percentile:
                raise ValueError(
                    f"{sample['label']}: expected Nx percentile N{expected_percentile}"
                )
            value = point.get("y")
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value <= 0):
                raise ValueError(
                    f"{sample['label']}: invalid read length at N{expected_percentile}"
                )
            points.append({"x": expected_percentile, "y": round(value, 3)})
        series.append({
            "name": sample["label"],
            "color": sample["color"],
            "data": points,
        })
    return {"title": "Read Length Nx", "samples": series}


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
# HTML TEMPLATES
# ---------------------------------------------------------------------------
# Templates are filled with report_common.fill(): %%NAME%% tokens are replaced verbatim, so
# CSS and JS braces need no escaping.

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>%%TITLE%%</title>
<meta name="description" content="Interactive ONT sequencing run comparison report.">
%%D3_SOURCE%%
<style>
  *,*::before,*::after{box-sizing:border-box;}
  body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
        max-width:1160px;margin:0 auto;padding:24px 28px;color:#2d3748;
        line-height:1.55;background:#f7f8fa;}
  h2{color:#1a202c;border-bottom:2px solid #e2e8f0;padding-bottom:8px;margin-top:40px;}
  h3{color:#2d3748;margin-top:28px;font-size:1.05rem;}
  table{border-collapse:collapse;width:100%;margin:0;font-size:13.5px;background:#fff;}
  th{background:#edf2f7;color:#2d3748;font-weight:600;text-align:left;
      padding:9px 10px;border-bottom:2px solid #cbd5e0;
      cursor:pointer;user-select:none;white-space:nowrap;}
  th:hover{background:#e2e8f0;}
  .sort-arrow{font-size:10px;color:#4a5568;margin-left:3px;}
  td{padding:8px 10px 14px;border-bottom:1px solid #e9ecef;position:relative;vertical-align:middle;
     white-space:nowrap;}
  tbody tr:last-child td{border-bottom:none;}
  tbody tr:hover td{background:rgba(59,130,246,.04);}
  .bar-track{position:absolute;bottom:4px;left:10px;right:10px;height:3px;
              background:rgba(59,130,246,.12);border-radius:99px;overflow:hidden;}
  .bar-fill{position:absolute;inset:0;height:100%;border-radius:99px;
             background:#3b82f6;
             transition:width .35s cubic-bezier(.4,0,.2,1);}
  tfoot tr td{background:#edf2f7 !important;font-weight:600;border-top:2px solid #cbd5e0;
               color:#2d3748;padding:8px 10px;white-space:nowrap;border-bottom:none;}
  tfoot td.mean{font-style:italic;font-weight:500;color:#4a5568;}
  /* Column groups above the column names */
  thead tr.groups th{cursor:default;background:#e2e8f0;color:#4a5568;font-size:11px;
                     text-transform:uppercase;letter-spacing:.04em;padding:5px 12px;
                     border-bottom:1px solid #cbd5e0;text-align:center;}
  thead tr.groups th + th{border-left:2px solid #f7f8fa;}
  th.group-start,td.group-start{border-left:2px solid #edf2f7;}
  th.nosort{cursor:default;}
  th.nosort:hover{background:#edf2f7;}
  /* Sample column stays visible while the table scrolls sideways */
  #summary-table th:first-child,#summary-table td:first-child{position:sticky;left:0;z-index:2;
     box-shadow:inset -1px 0 0 #e2e8f0;}
  #summary-table tbody td:first-child{background:#fff;}
  #summary-table thead tr.groups th:first-child{background:#e2e8f0;}
  /* QC threshold chips */
  .qc{display:inline-block;padding:1px 7px;border-radius:99px;font-variant-numeric:tabular-nums;}
  .qc-good{background:#dcfce7;color:#166534;}
  .qc-warn{background:#fef3c7;color:#92400e;}
  .qc-bad{background:#fee2e2;color:#991b1b;}
  td.spark{padding:4px 10px;line-height:0;}
  .sub{display:block;font-size:11px;color:#718096;line-height:1.2;}
  td.spark svg{display:block;}
  .toggle{font-size:12px;color:#4a5568;display:inline-flex;align-items:center;gap:5px;cursor:pointer;}
  .head-actions{display:flex;align-items:center;gap:14px;}
  /* Rank mode: each metric column coloured best -> worst; bars and chips off */
  .rank-mode .bar-track{display:none;}
  .rank-mode .qc{background:none;color:inherit;padding:0;}
  .table-wrap{overflow-x:auto;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.07);margin:10px 0 6px;}
  .section-head{display:flex;align-items:center;justify-content:space-between;gap:12px;
                flex-wrap:wrap;margin-top:28px;}
  .section-head h3{margin:0;}
  .btn{font:inherit;font-size:12px;padding:4px 12px;border:1px solid #cbd5e0;border-radius:6px;
       background:#fff;color:#2d3748;cursor:pointer;}
  .btn:hover{background:#edf2f7;}
  .table-note{font-size:12px;color:#718096;margin:0 0 8px;}
  .plot-container{margin:20px 0;background:#fff;border-radius:8px;
                   box-shadow:0 1px 4px rgba(0,0,0,.08);border:1px solid #e2e8f0;overflow:hidden;}
  .plot-container iframe{border:none;display:block;width:100%;}
  .footer{margin-top:48px;border-top:1px solid #e2e8f0;padding-top:14px;
           color:#a0aec0;font-size:12px;text-align:center;}
  @media (max-width:640px){
    body{padding:16px;}
    h2{margin-top:16px;font-size:1.3rem;}
  }
</style>
</head>
<body>

<h2>%%TITLE%%</h2>

<div class="section-head">
  <h3>Metrics Summary &amp; Pore Retention</h3>
  <div class="head-actions">
    <label class="toggle" data-tip="<div class='tt-title'>Colour by rank</div><div class='tt-text'>Shade each metric column from best (green) to worst (red) flow cell. Blocking: lower is better.</div>">
      <input type="checkbox" id="rank-toggle" onchange="applyRankColours()"> Colour by rank
    </label>
    <button class="btn" onclick="downloadTableCsv('summary-table', 'summary')">Download CSV</button>
  </div>
</div>
<div class="table-wrap">
<table id="summary-table">
  <thead>
%%SUMMARY_HEADER%%
  </thead>
  <tbody id="summary-tbody">
%%SUMMARY_ROWS%%
  </tbody>
  <tfoot>
%%TOTAL_ROW%%
  </tfoot>
</table>
</div>
<p class="table-note">%%TABLE_NOTE%%</p>
<script>
%%TABLE_CSV_JS%%

  // Sortable columns live in the bottom header row.
  const SUMMARY_THS = Array.from(document.querySelectorAll('#summary-table thead tr:last-child th'));

  function sortTable(th) {
    const tbody = th.closest('table').querySelector('tbody');
    const colIndex = SUMMARY_THS.indexOf(th);
    const type = th.dataset.type;
    const dir = th.dataset.dir === 'asc' ? 'desc' : 'asc';
    SUMMARY_THS.forEach(h => {
      h.dataset.dir = '';
      const arrow = h.querySelector('.sort-arrow');
      if (arrow) arrow.textContent = '';
    });
    th.dataset.dir = dir;
    th.querySelector('.sort-arrow').textContent = dir === 'asc' ? '\u25B2' : '\u25BC';
    const rows = Array.from(tbody.querySelectorAll('tr'));
    rows.sort((a, b) => {
      const av = a.children[colIndex].dataset.value;
      const bv = b.children[colIndex].dataset.value;
      let cmp;
      if (type === 'number') {
        const an = av === '' ? -Infinity : parseFloat(av);
        const bn = bv === '' ? -Infinity : parseFloat(bv);
        cmp = an - bn;
      } else {
        cmp = av.localeCompare(bv);
      }
      return dir === 'asc' ? cmp : -cmp;
    });
    rows.forEach(r => tbody.appendChild(r));
  }

  // RANK COLOURS: per column with data-better, red (worst) -> amber -> green
  // (best) across the flow cells. Columns where all values tie stay plain.
  const RANK_STOPS = [[254, 226, 226], [254, 243, 199], [220, 252, 231]];
  function rankColour(t) {
    const [a, b, u] = t < 0.5 ? [RANK_STOPS[0], RANK_STOPS[1], t * 2] : [RANK_STOPS[1], RANK_STOPS[2], t * 2 - 1];
    return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * u)).join(',')})`;
  }
  function applyRankColours() {
    const on = document.getElementById('rank-toggle').checked;
    const table = document.getElementById('summary-table');
    table.classList.toggle('rank-mode', on);
    const rows = Array.from(document.querySelectorAll('#summary-tbody tr'));
    SUMMARY_THS.forEach((th, ci) => {
      const better = th.dataset.better;
      if (!better) return;
      const vals = rows.map(r => parseFloat(r.children[ci].dataset.value));
      const min = Math.min(...vals), max = Math.max(...vals);
      rows.forEach((r, ri) => {
        const td = r.children[ci];
        if (!on || !(max > min) || isNaN(vals[ri])) { td.style.background = ''; return; }
        const t = (vals[ri] - min) / (max - min);
        td.style.background = rankColour(better === 'low' ? 1 - t : t);
      });
    });
  }

  // DATA BARS: thin inset pill under amounts that start at zero (data-bar)
  const BAR_COL_INDICES = SUMMARY_THS.map((th, i) => th.dataset.bar === 'true' ? i : -1).filter(i => i >= 0);

  function ensureTrack(td) {
    let track = td.querySelector('.bar-track');
    if (!track) {
      track = document.createElement('div');
      track.className = 'bar-track';
      const fill = document.createElement('div');
      fill.className = 'bar-fill';
      track.appendChild(fill);
      td.appendChild(track);
    }
    return track.querySelector('.bar-fill');
  }

  function applyDataBars() {
    const rows = Array.from(document.querySelectorAll('#summary-tbody tr'));
    BAR_COL_INDICES.forEach(ci => {
      let maxVal = 0;
      rows.forEach(r => {
        const v = parseFloat(r.children[ci]?.dataset.value);
        if (!isNaN(v) && v > maxVal) maxVal = v;
      });
      rows.forEach(r => {
        const td = r.children[ci];
        if (!td) return;
        const fill = ensureTrack(td);
        const v = parseFloat(td.dataset.value);
        fill.style.width = (!isNaN(v) && maxVal > 0)
          ? `${Math.round((v / maxVal) * 100)}%`
          : '0%';
      });
    });
  }

  applyDataBars();
</script>

%%NX_SECTION%%

<h3>Active Pores Over Time (Pore Stability)</h3>
<div class="plot-container">
  <iframe srcdoc="%%ACTIVE_PORES_IFRAME%%" height="720" frameborder="0" scrolling="no"></iframe>
</div>

<h3>Cumulative Yield Over Time</h3>
<div class="plot-container">
  <iframe srcdoc="%%YIELD_IFRAME%%" height="520" frameborder="0" scrolling="no"></iframe>
</div>

<h3>Average Pore Activity</h3>
<div class="plot-container">
  <iframe srcdoc="%%PORE_ACTIVITY_IFRAME%%" height="820" frameborder="0" scrolling="no"></iframe>
</div>

%%FOOTER%%
%%PAGE_TIP%%
</body>
</html>
"""


LINE_CHART_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>%%TITLE%%</title>
%%CHART_HEAD%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title">%%TITLE%%</h2>
%%EXPORT_BUTTONS%%
    <div class="toggle-bar">
      <button class="toggle-btn active" data-mode="yield" onclick="setMode('yield')">Yield Over Time</button>
      <button class="toggle-btn" data-mode="rate" onclick="setMode('rate')">Rate (Gb/min)</button>
      <button class="toggle-btn" data-mode="perpore" onclick="setMode('perpore')">Rate per Pore (b/s)</button>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = %%PLOT_DATA%%;
    const CHART_HEIGHT = %%CHART_HEIGHT%%;
    const Y_LABEL = %%Y_LABEL%%;
    unpack(plotData.samples, ['data', 'estimated', 'pores']);

    // First derivative of cumulative yield: windowed rate in Gb/min.
    // Linearly interpolates the cumulative curve at regular BIN_MIN intervals,
    // then computes slope = ΔY / BIN_MIN — giving a smooth, readable signal
    // instead of noisy point-by-point differences on dense throughput data.
    const BIN_MIN = 100;
    function interpY(data, t) {
      if (t <= data[0].x) return data[0].y;
      if (t >= data[data.length-1].x) return data[data.length-1].y;
      let lo = 0, hi = data.length - 1;
      while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (data[mid].x <= t) lo = mid; else hi = mid;
      }
      const frac = (t - data[lo].x) / (data[hi].x - data[lo].x);
      return data[lo].y + frac * (data[hi].y - data[lo].y);
    }
    // Mean number of sequencing-capable pores (strand + adapter + pore) per
    // minute in (t0, t1] — the denominator for the per-pore rate.
    function meanPores(pores, t0, t1) {
      let sum = 0, n = 0;
      for (const p of pores) {
        if (p.x > t0 && p.x <= t1) { sum += p.y; n++; }
      }
      return n ? sum / n : 0;
    }
    // Windows with fewer pores than this fraction of the run's peak are
    // dropped from the per-pore rate: the ratio is just noise near zero.
    const MIN_PORE_FRAC = 0.05;
    plotData.samples.forEach(s => {
      if (s.data.length < 2) { s.rateData = []; s.perPoreData = []; return; }
      const xMax = s.data[s.data.length - 1].x;
      const pores = s.pores || [];
      const minPores = MIN_PORE_FRAC * d3.max(pores, p => p.y) || 1;
      const pts = [], perPore = [];
      for (let t = BIN_MIN; t <= xMax + BIN_MIN / 2; t += BIN_MIN) {
        const t1 = Math.min(t, xMax);
        const t0 = t1 - BIN_MIN;
        const rate = Math.max(0, (interpY(s.data, t1) - interpY(s.data, t0)) / BIN_MIN);
        pts.push({x: t1, y: rate});
        const n = meanPores(pores, t0, t1);
        if (n >= minPores) {
          const est = s.estimated || s.data;
          const estRate = Math.max(0, (interpY(est, t1) - interpY(est, t0)) / BIN_MIN);
          perPore.push({x: t1, y: Math.round(estRate * 1e9 / 60 / n)});
        }
      }
      s.rateData = pts;
      s.perPoreData = perPore;
    });

    const MODES = {
      yield:   {key: 'data', label: Y_LABEL, title: plotData.title,
                fmt: v => v.toFixed(2) + ' Gb'},
      rate:    {key: 'rateData', label: 'Rate (Gb/min)', title: 'Yield Rate Over Time (Gb/min)',
                fmt: v => v.toFixed(3) + ' Gb/min'},
      perpore: {key: 'perPoreData', label: 'Bases / s / sequencing pore',
                title: 'Estimated Yield Rate per Sequencing Pore (strand + adapter + pore)',
                fmt: v => v.toLocaleString() + ' b/s'},
    };
    let mode = 'yield'; // 'yield' | 'rate' | 'perpore'
    let draw = () => {};

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const totalW = chartWidth();
      const legend = legendLayout(totalW, plotData.samples.length);
      const margin = {top: 20, right: legend.right, bottom: 50, left: 65};
      const W = totalW - margin.left - margin.right;
      const H = CHART_HEIGHT - margin.top - margin.bottom;
      const svg = appendSvg(totalW, margin, H, legend, 56);

      const allPoints = plotData.samples.flatMap(s => s.data);
      const x = d3.scaleLinear().domain(d3.extent(allPoints, d => d.x / 60)).range([0, W]);
      const y = d3.scaleLinear().range([H, 0]);

      const gGrid  = svg.append('g').attr('class', 'grid');
      const gXAxis = svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`);
      const gYAxis = svg.append('g').attr('class', 'axis');
      gXAxis.call(hoursAxis(x, W));
      axisLabel(svg, 'Experiment Time (hours)', {x: W / 2, y: H + 42});
      const yLabel = axisLabel(svg, '', {transform: 'rotate(-90)', x: -H / 2, y: -52});

      const lineGen = d3.line().x(d => x(d.x / 60)).curve(d3.curveMonotoneX);

      // One path + dots group per sample (updated on redraw)
      plotData.samples.forEach((s, i) => {
        const g = svg.append('g').attr('class', `sample-group sample-${i}`);
        g.append('path').attr('class', 'line').attr('stroke', s.color);
        g.append('g').attr('class', 'dots');
      });

      const hoverLine = svg.append('line').attr('class', 'hover-line').attr('y1', 0).attr('y2', H);
      const overlay = svg.append('rect')
          .attr('width', W).attr('height', H).attr('fill', 'none').attr('pointer-events', 'all');

      drawLegend(svg, plotData.samples, W, H, legend, 56);

      draw = function () {
        const m = MODES[mode];
        document.getElementById('chart-title').textContent = m.title;

        const activePoints = plotData.samples.flatMap(s => s[m.key]);
        y.domain([0, (d3.max(activePoints, d => d.y) || 1) * 1.05]);
        gGrid.call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
        gYAxis.call(d3.axisLeft(y));
        yLabel.text(m.label);
        lineGen.y(d => y(d.y));

        plotData.samples.forEach((s, i) => {
          const g = svg.select(`.sample-${i}`);
          g.select('path').datum(s[m.key]).attr('d', lineGen);
          // Dot markers — shown only in the rate modes
          g.select('.dots').selectAll('circle')
             .data(mode === 'yield' ? [] : s[m.key])
             .join('circle')
             .attr('r', 3.5)
             .attr('cx', d => x(d.x / 60))
             .attr('cy', d => y(d.y))
             .attr('fill', s.color)
             .attr('stroke', '#fff')
             .attr('stroke-width', 1.5);
        });

        overlay
          .on('mousemove', function (event) {
            const [mx] = d3.pointer(event);
            const t = x.invert(mx);
            hoverLine.attr('x1', mx).attr('x2', mx).style('opacity', 1);
            const rows = [];
            plotData.samples.forEach(s => {
              if (s.hidden) return;
              const data = s[m.key];
              const pt = data[d3.bisector(d => d.x).left(data, t * 60)] || data[data.length - 1];
              if (pt) rows.push([s.color, s.name, pt.y]);
            });
            rows.sort((a, b) => b[2] - a[2]);
            showTooltip(event, tipHtml(`${t.toFixed(1)} h`, rows.map(([c, n, v]) => [c, n, m.fmt(v)])));
          })
          .on('mouseout', () => { hoverLine.style('opacity', 0); hideTooltip(); });
      };

      draw();
      fitFrame();
    }

    function setMode(m) {
      mode = m;
      document.querySelectorAll('.toggle-btn').forEach(b =>
        b.classList.toggle('active', b.dataset.mode === m));
      draw();
    }

    build();
    onResize(build);
  </script>
</body>
</html>"""


NX_CHART_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>%%TITLE%%</title>
%%CHART_HEAD%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title">%%TITLE%%</h2>
%%EXPORT_BUTTONS%%
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = %%PLOT_DATA%%;
    const CHART_HEIGHT = %%CHART_HEIGHT%%;
    const percentiles = Array.from({length: 19}, (_, i) => (i + 1) * 5);

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const totalW = chartWidth();
      const legend = legendLayout(totalW, plotData.samples.length);
      const margin = {top: 20, right: legend.right, bottom: 56, left: 70};
      const W = totalW - margin.left - margin.right;
      const H = CHART_HEIGHT - margin.top - margin.bottom;
      const svg = appendSvg(totalW, margin, H, legend, 62);

      const x = d3.scaleLinear().domain([5, 95]).range([0, W]);
      const yMax = d3.max(plotData.samples.flatMap(s => s.data), d => d.y) || 1;
      const y = d3.scaleLinear().domain([0, yMax * 1.05]).nice().range([H, 0]);
      // Label every other percentile when space is tight.
      const ticks = W < 500 ? percentiles.filter(p => p % 10 === 0) : percentiles;
      const xAxis = svg.append('g').attr('class', 'axis')
          .attr('transform', `translate(0,${H})`)
          .call(d3.axisBottom(x).tickValues(ticks).tickFormat(d => 'N' + d));
      svg.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickFormat(d => d + ' kb'));
      svg.append('g').attr('class', 'grid').call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
      xAxis.selectAll('text').attr('transform', 'rotate(-35)')
          .style('text-anchor', 'end').attr('dx', '-0.4em').attr('dy', '0.15em');
      axisLabel(svg, 'Nx percentile', {x: W / 2, y: H + 48});
      axisLabel(svg, 'Read length (kb)', {transform: 'rotate(-90)', x: -H / 2, y: -55});

      const line = d3.line().x(d => x(d.x)).y(d => y(d.y)).curve(d3.curveMonotoneX);
      plotData.samples.forEach((sample, i) => {
        const group = svg.append('g').attr('class', `sample-group sample-${i}`);
        group.append('path').datum(sample.data).attr('class', 'line')
            .attr('stroke', sample.color).attr('d', line);
        group.selectAll('circle').data(sample.data).join('circle')
            .attr('cx', d => x(d.x)).attr('cy', d => y(d.y)).attr('r', 3)
            .attr('fill', sample.color).attr('stroke', '#fff').attr('stroke-width', 1);
      });

      const hoverLine = svg.append('line').attr('class', 'hover-line').attr('y1', 0).attr('y2', H);
      svg.append('rect').attr('width', W).attr('height', H).attr('fill', 'none')
          .attr('pointer-events', 'all')
          .on('mousemove', function (event) {
            const [mx] = d3.pointer(event);
            const nearest = percentiles.reduce((best, p) =>
              Math.abs(x(p) - mx) < Math.abs(x(best) - mx) ? p : best, percentiles[0]);
            const isN50 = nearest === 50;
            hoverLine.attr('x1', x(nearest)).attr('x2', x(nearest))
                .style('stroke', isN50 ? '#dc2626' : '#888')
                .style('stroke-width', isN50 ? '3px' : '1px')
                .style('stroke-dasharray', isN50 ? '6,3' : '4,4')
                .style('opacity', 1);
            const rows = plotData.samples
              .filter(sample => !sample.hidden)
              .map(sample => [sample.color, sample.name, (sample.data.find(d => d.x === nearest) || {}).y])
              .filter(r => r[2] !== undefined)
              .sort((a, b) => b[2] - a[2])
              .map(([c, n, v]) => [c, n, `${v.toFixed(1)} kb`]);
            showTooltip(event, tipHtml(`N${nearest}`, rows,
              `${nearest}% of all bases are in reads at least this long`));
          })
          .on('mouseout', () => { hoverLine.style('opacity', 0); hideTooltip(); });

      drawLegend(svg, plotData.samples, W, H, legend, 62);
      fitFrame();
    }

    build();
    onResize(build);
  </script>
</body>
</html>"""


ACTIVE_PORES_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>%%TITLE%%</title>
%%CHART_HEAD%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title">%%TITLE%%</h2>
%%EXPORT_BUTTONS%%
    <div class="toggle-bar">
      <button class="toggle-btn active" data-mode="abs" onclick="setMode('abs')">Absolute</button>
      <button class="toggle-btn" data-mode="rel" onclick="setMode('rel')">Relative (%)</button>
      <button class="toggle-btn" data-mode="occ" onclick="setMode('occ')">Occupancy (%)</button>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = %%PLOT_DATA%%;
    const CHART_HEIGHT = %%CHART_HEIGHT%%;

    unpack(plotData.samples, ['data', 'occupancy']);

    // 100% = peak pore count for the relative view
    plotData.samples.forEach(s => {
      const peak = d3.max(s.data, d => d.y) || 1;
      s.relData = s.data.map(d => ({x: d.x, y: +(d.y / peak * 100).toFixed(2)}));
      s.occData = s.occupancy || [];
    });

    const pct = v => v.toFixed(1) + '%';
    const MODES = {
      abs: {key: 'data', label: 'Active Pores', title: plotData.title,
            fmt: v => v.toLocaleString(), tick: d => d.toLocaleString()},
      rel: {key: 'relData', label: 'Active Pores (% of peak)', title: 'Active Pores Over Time (% of peak)',
            fmt: pct, tick: d => d + '%'},
      occ: {key: 'occData', label: 'Pore occupancy (%)',
            title: 'Pore Occupancy: (strand + adapter) / (strand + adapter + pore)',
            fmt: pct, tick: d => d + '%', yMax: 100},
    };
    let mode = 'abs'; // 'abs' | 'rel' | 'occ'
    let render = () => {};

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const totalW = chartWidth();
      const legend = legendLayout(totalW, plotData.samples.length);
      const margin = {top: 20, right: legend.right, bottom: 50, left: 65};
      const W = totalW - margin.left - margin.right;
      const H = CHART_HEIGHT - margin.top - margin.bottom;
      const svg = appendSvg(totalW, margin, H, legend, 56);

      const allPoints = plotData.samples.flatMap(s => s.data);
      const x = d3.scaleLinear().domain(d3.extent(allPoints, d => d.x / 60)).range([0, W]);
      const y = d3.scaleLinear().range([H, 0]);

      const gridG  = svg.append('g').attr('class', 'grid');
      const xAxisG = svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`);
      const yAxisG = svg.append('g').attr('class', 'axis');
      xAxisG.call(hoursAxis(x, W));
      axisLabel(svg, 'Experiment Time (hours)', {x: W / 2, y: H + 42});
      const yLabel = axisLabel(svg, '', {transform: 'rotate(-90)', x: -H / 2, y: -52});

      const lineGen = d3.line().x(d => x(d.x / 60)).y(d => y(d.y)).curve(d3.curveMonotoneX);
      plotData.samples.forEach((s, i) => {
        svg.append('g').attr('class', `sample-group sample-${i}`)
           .append('path').attr('class', 'line').attr('stroke', s.color);
      });

      const hoverLine = svg.append('line').attr('class', 'hover-line').attr('y1', 0).attr('y2', H);
      svg.append('rect')
          .attr('width', W).attr('height', H).attr('fill', 'none').attr('pointer-events', 'all')
          .on('mousemove', function (event) {
            const [mx] = d3.pointer(event);
            const t = x.invert(mx);
            const m = MODES[mode];
            hoverLine.attr('x1', mx).attr('x2', mx).style('opacity', 1);
            const rows = [];
            plotData.samples.forEach(s => {
              if (s.hidden) return;
              const pts = s[m.key];
              const d = pts[d3.bisector(p => p.x).left(pts, t * 60)] || pts[pts.length - 1];
              if (d) rows.push([s.color, s.name, d.y]);
            });
            rows.sort((a, b) => b[2] - a[2]);
            showTooltip(event, tipHtml(`${t.toFixed(1)} h`, rows.map(([c, n, v]) => [c, n, m.fmt(v)])));
          })
          .on('mouseout', () => { hoverLine.style('opacity', 0); hideTooltip(); });

      drawLegend(svg, plotData.samples, W, H, legend, 56);

      render = function (animate) {
        const m = MODES[mode];
        document.getElementById('chart-title').textContent = m.title;
        const pts = plotData.samples.flatMap(s => s[m.key]);
        y.domain([0, m.yMax || (d3.max(pts, d => d.y) || 1) * 1.05]);
        const ms = animate ? 400 : 0;
        yAxisG.transition().duration(ms).call(d3.axisLeft(y).tickFormat(m.tick));
        gridG.transition().duration(ms).call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
        yLabel.text(m.label);
        plotData.samples.forEach((s, i) => {
          svg.select(`.sample-${i} path`).datum(s[m.key])
             .transition().duration(ms).attr('d', lineGen);
        });
      };

      render(false);
      fitFrame();
    }

    function setMode(m) {
      mode = m;
      document.querySelectorAll('.toggle-btn').forEach(b =>
        b.classList.toggle('active', b.dataset.mode === m));
      render(true);
    }

    build();
    onResize(build);
  </script>
</body>
</html>"""


STACKED_BAR_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>Average Pore Activity</title>
%%CHART_HEAD%%
  <style>
    .grid line{stroke:#eee;stroke-dasharray:none;}
    p.subtitle{text-align:center;color:#666;font-size:13px;margin:0 0 16px;}
    .legend-container{display:flex;flex-wrap:wrap;justify-content:center;gap:8px;
                       margin-top:18px;padding:12px;background:#fcfcfc;border-radius:8px;
                       border:1px solid #eee;}
    .legend-container .legend-item{display:flex;align-items:center;font-size:12px;color:#444;
                  padding:5px 10px;border-radius:5px;transition:all .2s;user-select:none;
                  border:1px solid transparent;}
    .legend-container .legend-item:hover{background:#f0f0f0;border-color:#ddd;}
    .legend-item.hidden{opacity:.35;}
    .legend-item.hidden .legend-text{text-decoration:line-through;color:#888;}
    .legend-color{width:13px;height:13px;margin-right:7px;border-radius:3px;
                   border:1px solid rgba(0,0,0,.12);}
  </style>
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title">Average Pore Activity per Sample</h2>
%%EXPORT_BUTTONS%%
    <p class="subtitle">Click legend items to toggle states</p>
    <div id="chart"></div>
    <div class="legend-container" id="legend"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const plotData = %%PLOT_DATA%%;
    let visibleStates = [...plotData.states];
    let updateChart = () => {};

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const totalW = chartWidth();
      const margin = {top: 20, right: 20, bottom: 100, left: 65};
      const W = totalW - margin.left - margin.right;
      const H = 500 - margin.top - margin.bottom;

      const svg = d3.select('#chart').append('svg')
          .attr('width',  totalW)
          .attr('height', H + margin.top + margin.bottom)
        .append('g')
          .attr('transform', `translate(${margin.left},${margin.top})`);

      const x = d3.scaleBand().domain(plotData.samples.map(d => d.Sample)).range([0, W]).padding(0.3);
      const y = d3.scaleLinear().range([H, 0]);
      const stack = d3.stack();

      svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x)).selectAll('text')
         .attr('transform', 'rotate(-30)').style('text-anchor', 'end').attr('dy', '0.8em');

      const gridGroup = svg.append('g').attr('class', 'grid');
      const yAxisGroup = svg.append('g').attr('class', 'axis');
      axisLabel(svg, '% Time in State', {transform: 'rotate(-90)', x: -H / 2, y: -55});

      updateChart = function () {
        const totals = plotData.samples.map(d => visibleStates.reduce((s, k) => s + (d[k] || 0), 0));
        const maxTotal = d3.max(totals) || 100;
        y.domain([0, Math.max(100, maxTotal * 1.05)]);
        yAxisGroup.transition().duration(400).call(d3.axisLeft(y).tickFormat(d => d + '%'));
        gridGroup.transition().duration(400).call(d3.axisLeft(y).tickSize(-W).tickFormat(''));

        stack.keys(visibleStates);
        const layers = svg.selectAll('.layer').data(stack(plotData.samples), d => d.key);
        layers.exit().remove();
        const allLayers = layers.enter().append('g').attr('class', 'layer')
            .attr('fill', d => plotData.colors[d.key] || '#ccc')
            .merge(layers);

        let barTip = '';
        const rects = allLayers.selectAll('rect').data(d => d, d => d.data.Sample);
        rects.exit().remove();
        rects.enter().append('rect')
             .attr('x', d => x(d.data.Sample)).attr('width', x.bandwidth())
             .attr('y', H).attr('height', 0)
             .merge(rects)
             .on('mouseover', function (event) {
               const key = d3.select(this.parentNode).datum().key;
               const color = plotData.colors[key] || '#ccc';
               barTip = tipHtml(`<span class="tt-sw" style="background:${color}"></span>${escTip(key)}`,
                 plotData.samples.map(s => [null, s.Sample, `${(s[key] || 0).toFixed(2)}%`]),
                 '% of time in this state');
               showTooltip(event, barTip);
               d3.select(this).attr('opacity', 0.7);
             })
             .on('mousemove', event => showTooltip(event, barTip))
             .on('mouseout', function () {
               hideTooltip();
               d3.select(this).attr('opacity', 1);
             })
             .transition().duration(400)
             .attr('y', d => y(d[1])).attr('height', d => Math.max(0, y(d[0]) - y(d[1])));
      };

      updateChart();
      fitFrame();
    }

    // Export: the state legend is HTML, so redraw the visible states into the SVG.
    function stateLegend(svg, w, h) {
      const cols = Math.max(1, Math.floor((w - 40) / 210));
      const g = document.createElementNS(SVG_NS, 'g');
      g.setAttribute('transform', `translate(40,${h})`);
      visibleStates.forEach((state, i) => {
        const cx = (i % cols) * 210, cy = Math.floor(i / cols) * 20;
        const box = document.createElementNS(SVG_NS, 'rect');
        box.setAttribute('x', cx);
        box.setAttribute('y', cy);
        box.setAttribute('width', 12);
        box.setAttribute('height', 12);
        box.setAttribute('fill', plotData.colors[state] || '#ccc');
        const label = document.createElementNS(SVG_NS, 'text');
        label.setAttribute('x', cx + 18);
        label.setAttribute('y', cy + 10);
        label.setAttribute('style', 'font-size:12px;fill:#444;font-family:'
          + getComputedStyle(document.body).fontFamily);
        label.textContent = state;
        g.append(box, label);
      });
      svg.appendChild(g);
      return Math.ceil(visibleStates.length / cols) * 20 + 10;
    }

    const legendContainer = d3.select('#legend');
    plotData.states.forEach(state => {
      const item = legendContainer.append('div').attr('class', 'legend-item')
          .on('click', function () {
            const idx = visibleStates.indexOf(state);
            if (idx > -1) {
              if (visibleStates.length > 1) {
                visibleStates.splice(idx, 1);
                d3.select(this).classed('hidden', true);
              }
            } else {
              visibleStates.push(state);
              visibleStates.sort((a, b) => plotData.states.indexOf(a) - plotData.states.indexOf(b));
              d3.select(this).classed('hidden', false);
            }
            updateChart();
          });
      item.append('div').attr('class', 'legend-color').style('background-color', plotData.colors[state] || '#ccc');
      item.append('span').attr('class', 'legend-text').text(state);
    });

    build();
    onResize(build);
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# RENDERING HELPERS
# ---------------------------------------------------------------------------

def _render_chart_iframe(template, plot_data, chart_height=0, decorate="", **extra):
    return render_chart_iframe(template, plot_data, plot_data.get("title", ""), decorate,
                               chart_height=chart_height, **extra)


def _pack_series(plot_data, keys):
    """Copy of plot_data with {x, y} point lists packed as [x, y] pairs.

    Inside a srcdoc attribute every quote becomes &quot;, so dense time series
    shrink several-fold without their "x"/"y" keys. Charts unpack with unpack().
    """
    samples = [{k: [[p["x"], p["y"]] for p in v] if k in keys else v
                for k, v in s.items()}
               for s in plot_data["samples"]]
    return {**plot_data, "samples": samples}


def _render_line_chart_iframe(plot_data, chart_height, y_label):
    packed = _pack_series(plot_data, {"data", "estimated", "pores"})
    return _render_chart_iframe(LINE_CHART_IFRAME, packed, chart_height,
                                y_label=json.dumps(y_label))


def _render_nx_iframe(plot_data, chart_height):
    return _render_chart_iframe(NX_CHART_IFRAME, plot_data, chart_height)


def _render_active_pores_iframe(plot_data, chart_height):
    packed = _pack_series(plot_data, {"data", "occupancy"})
    return _render_chart_iframe(ACTIVE_PORES_IFRAME, packed, chart_height)


def _render_stacked_bar_iframe(plot_data):
    return _render_chart_iframe(STACKED_BAR_IFRAME, plot_data, decorate=", stateLegend")


# QC thresholds for colouring summary-table cells: green when at least `good`,
# amber when at least `warn`, else red ("low" metrics: at most, i.e. reversed).
# Starting points — adjust to your lab's expectations.
QC_THRESHOLDS = {
    "pass_rate":     {"good": 85, "warn": 75, "better": "high"},
    "occupancy_pct": {"good": 85, "warn": 70, "better": "high"},
    "blocking_pct":  {"good": 5,  "warn": 10, "better": "low"},
}

# Summary table columns, grouped under headings.
#   kind: text | number | spark;  bar: zero-based amount, gets a data bar;
#   qc: coloured against QC_THRESHOLDS;  better: direction for rank colours;
#   total: sum | pooled | mean | label (how the Total row is filled)
SUMMARY_GROUPS = [
    ("Flow cell", [
        # The flow cell ID is shown under the label when it differs; CSV gets both.
        {"key": "label", "label": "Sample", "kind": "text", "total": "label",
         "csv_cols": "Sample|Flow Cell ID"},
    ]),
    ("Output", [
        {"key": "yield_gb", "label": "Yield (Gb)", "kind": "number", "fmt": "{:,.2f}",
         "bar": True, "better": "high", "total": "sum"},
        {"key": "passed_reads", "label": "Passed Reads", "kind": "number", "fmt": "{:,.0f}",
         "bar": True, "better": "high", "total": "sum"},
    ]),
    ("Quality", [
        {"key": "pass_rate", "label": "Pass Rate", "csv_cols": "Pass Rate (%)", "kind": "number", "fmt": "{:.2f}%",
         "qc": True, "better": "high", "total": "pooled"},
    ]),
    ("Pore health", [
        {"key": "avg_active", "label": "Avg Pores", "csv_cols": "Avg Active Pores", "kind": "number",
         "fmt": "{:,.0f}", "bar": True, "better": "high", "total": "mean",
         "title": "Run-average number of active pores"},
        {"key": "yield_per_pore_mb", "label": "Mb / Pore", "csv_cols": "Yield / Pore (Mb)", "kind": "number", "fmt": "{:.2f}",
         "bar": True, "better": "high", "total": "mean",
         "title": "Yield divided by the run-average number of sequencing-capable pores (strand + adapter + pore)"},
        {"key": "occupancy_pct", "label": "Occupancy", "csv_cols": "Occupancy (%)", "kind": "number", "fmt": "{:.2f}%",
         "qc": True, "better": "high", "total": "mean",
         "title": "Of the pores able to sequence, the share of time spent in strand or adapter: "
                  "(strand + adapter) / (strand + adapter + pore), time-weighted over the run"},
        {"key": "blocking_pct", "label": "Blocking", "csv_cols": "Blocking (%)", "kind": "number", "fmt": "{:.2f}%",
         "qc": True, "better": "low", "total": "mean",
         "title": "Share of active pore time spent blocked (multiple or saturated)"},
    ]),
    ("Retention", [
        {"key": "retention_curve", "label": "Active pores, % of peak", "kind": "spark",
         "title": "Active pores as % of their peak over the run (60-min rolling median); "
                  "dashed line = 50%, dot = half-life. Time axis shared across rows."},
        {"key": "half_life_h", "label": "Half-life (h)", "kind": "number", "fmt": "{:.1f}",
         "better": "high", "total": "mean",
         "title": "Hours until active pores first fall below 50% of their peak; "
                  "\u2265 when they never did before the run ended"},
    ]),
]
SUMMARY_COLUMNS = [dict(col, group_start=(i == 0))
                   for _, cols in SUMMARY_GROUPS for i, col in enumerate(cols)]


def _qc_class(key, v):
    t = QC_THRESHOLDS[key]
    ok = (lambda lim: v >= lim) if t["better"] == "high" else (lambda lim: v <= lim)
    return "qc-good" if ok(t["good"]) else "qc-warn" if ok(t["warn"]) else "qc-bad"


def _sparkline(curve, max_minute, color, half_life_h, censored, w=130, h=34):
    """Inline SVG retention curve (% of peak) on a time axis shared by all rows."""
    pad = 3
    def px(m): return pad + (w - 2 * pad) * m / max(max_minute, 1)
    def py(pct): return pad + (h - 2 * pad) * (1 - min(max(pct, 0), 100) / 100)
    points = " ".join(f"{px(m):.1f},{py(p):.1f}" for m, p in curve)
    marker = "" if censored else (
        f'<circle cx="{px(half_life_h * 60):.1f}" cy="{py(50):.1f}" r="2.6" fill="{color}" '
        f'stroke="#fff" stroke-width="1"/>')
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" '
            f'aria-label="Active pores retention curve">'
            f'<line x1="{pad}" x2="{w - pad}" y1="{py(50):.1f}" y2="{py(50):.1f}" '
            f'stroke="#cbd5e0" stroke-dasharray="2,2"/>'
            f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="1.5" '
            f'stroke-linejoin="round"/>{marker}</svg>')


def _th_attrs(col):
    attrs = []
    classes = ["group-start"] if col["group_start"] else []
    if col["kind"] == "spark":
        classes.append("nosort")
        attrs.append('data-csv-skip=""')
    else:
        attrs += [f'data-type="{col["kind"]}"', 'onclick="sortTable(this)"']
    if col.get("csv_cols"):
        attrs.append(f'data-csv-cols="{col["csv_cols"]}"')
    if col["kind"] == "number":
        attrs.append('style="text-align:right"')
    if col.get("bar"):
        attrs.append('data-bar="true"')
    if col.get("better"):
        attrs.append(f'data-better="{col["better"]}"')
    if col.get("title"):
        attrs.append(tip_attr(tip_html(col.get("csv_cols", col["label"]), text=col["title"])))
    if classes:
        attrs.append(f'class="{" ".join(classes)}"')
    return " ".join(attrs)


def _render_summary_header():
    groups = "".join(f'<th colspan="{len(cols)}">{name}</th>' for name, cols in SUMMARY_GROUPS)
    cols = "".join(
        f'<th {_th_attrs(c)}>{c["label"]}'
        + ('' if c["kind"] == "spark" else '<span class="sort-arrow"></span>') + '</th>'
        for c in SUMMARY_COLUMNS)
    return f'    <tr class="groups">{groups}</tr>\n    <tr>{cols}</tr>'


QC_FLAG_NAMES = {"qc-good": "Good", "qc-warn": "Check", "qc-bad": "Poor"}


def _qc_tip(col, v):
    t = QC_THRESHOLDS[col["key"]]
    op = "\u2265" if t["better"] == "high" else "\u2264"
    name = col.get("csv_cols", col["label"])
    return tip_html(f"{name}: {col['fmt'].format(v)}", [
        (None, "Flag", QC_FLAG_NAMES[_qc_class(col["key"], v)]),
        ("#22c55e", "Good", f"{op} {t['good']}%"),
        ("#f59e0b", "Check", f"{op} {t['warn']}%"),
        ("#ef4444", "Poor", "otherwise"),
    ])


def _number_cell(col, v, extra_class="", tip="", text=None):
    classes = " ".join(c for c in ("group-start" if col["group_start"] else "", extra_class) if c)
    cls = f' class="{classes}"' if classes else ""
    text = col["fmt"].format(v) if text is None else text
    if col.get("qc"):
        text = f'<span class="qc {_qc_class(col["key"], v)}">{text}</span>'
        tip = tip or _qc_tip(col, v)
    tip_a = f" {tip_attr(tip)}" if tip else ""
    return f"<td{cls} style='text-align:right' data-value=\"{v}\"{tip_a}>{text}</td>"


def _retention_tip(s):
    curve = s["retention_curve"]
    if not curve:
        return tip_html(s["label"], text="No pore data")
    low_m, low_p = min(curve, key=lambda c: c[1])
    end_m, end_p = curve[-1]
    half = (f"\u2265 {s['half_life_h']:.1f} h (not reached)" if s["half_life_censored"]
            else f"{s['half_life_h']:.1f} h")
    rows = [(None, "Half-life (50% of peak)", half)]
    # A dip before the end (e.g. before a flush/reload) is worth its own row.
    if low_m != end_m:
        rows.append((None, "Lowest", f"{low_p:.0f}% at {low_m / 60:.0f} h"))
    rows.append((None, "Run end", f"{end_p:.0f}% at {end_m / 60:.0f} h"))
    return tip_html(f'<span class="tt-sw" style="background:{s["color"]}"></span>'
                    f'{html_lib.escape(s["label"])}', rows,
                    "Active pores as % of their peak (60-min rolling median)")


def _render_summary_rows(samples):
    max_minute = max((s["retention_curve"][-1][0] for s in samples if s["retention_curve"]), default=1)
    rows = []
    for s in samples:
        tds = []
        for col in SUMMARY_COLUMNS:
            gs = ' class="group-start"' if col["group_start"] else ""
            if col["key"] == "label":
                label, fc = s["label"], s.get("flow_cell_id") or ""
                sub = f'<span class="sub">{html_lib.escape(fc)}</span>' if fc and fc != label else ""
                tds.append(f'<td data-value="{html_lib.escape(label.lower())}" '
                           f'data-csv="{html_lib.escape(label)}|{html_lib.escape(fc)}">'
                           f'{html_lib.escape(label)}{sub}</td>')
            elif col["kind"] == "spark":
                tds.append(f'<td class="spark{" group-start" if col["group_start"] else ""}" '
                           f'{tip_attr(_retention_tip(s))}>'
                           + _sparkline(s["retention_curve"], max_minute, s["color"],
                                        s["half_life_h"], s["half_life_censored"]) + '</td>')
            elif col["key"] == "half_life_h" and s["half_life_censored"]:
                tds.append(_number_cell(col, s["half_life_h"], text=f'\u2265 {s["half_life_h"]:.1f}',
                                        tip=tip_html("Half-life not reached", text=(
                                            "Active pores never fell below 50% of their peak; "
                                            f"the run ended at {s['half_life_h']:.1f} h."))))
            else:
                tds.append(_number_cell(col, s[col["key"]]))
        rows.append("    <tr>" + "".join(tds) + "</tr>")
    return "\n".join(rows)


def _render_total_row(samples):
    n = len(samples)
    total_reads = sum(s["total_reads"] for s in samples)
    passed_reads = sum(s["passed_reads"] for s in samples)
    tds = []
    for col in SUMMARY_COLUMNS:
        gs = ' class="group-start"' if col["group_start"] else ""
        total = col.get("total")
        if total == "label":
            cells = f"{n} flow cell{'s' if n != 1 else ''}"
            tds.append(f'<td data-csv="Total|{cells}"><strong>Total</strong>'
                       f'<span class="sub">{cells}</span></td>')
        elif total == "sum":
            v = round(sum(s[col["key"]] for s in samples), 2)
            tds.append(f"<td{gs} style='text-align:right' data-value=\"{v}\">{col['fmt'].format(v)}</td>")
        elif total == "pooled":
            v = round(100 * passed_reads / max(total_reads, 1), 2)
            tds.append(f"<td{gs} style='text-align:right' data-value=\"{v}\">{col['fmt'].format(v)}</td>")
        elif total == "mean":
            v = round(sum(s[col["key"]] for s in samples) / max(n, 1), 2)
            text = col["fmt"].format(v)
            if col["key"] == "half_life_h" and any(s["half_life_censored"] for s in samples):
                text = "\u2265 " + text
            cls = "mean group-start" if col["group_start"] else "mean"
            tip = tip_attr(tip_html(f"Mean: {text}", text=f"Mean of {n} flow cells"))
            tds.append(f"<td class='{cls}' style='text-align:right' data-value=\"{v}\" "
                       f"{tip}>{text}</td>")
        else:
            tds.append(f"<td{gs}></td>")
    return "    <tr>" + "".join(tds) + "</tr>"


def _render_table_note():
    def rule(label, key, unit="%"):
        t = QC_THRESHOLDS[key]
        op = "\u2265" if t["better"] == "high" else "\u2264"
        return f'{label} {op}&nbsp;{t["good"]}{unit} / {op}&nbsp;{t["warn"]}{unit}'
    return (
        "Coloured values are QC flags: "
        '<span class="qc qc-good">good</span> <span class="qc qc-warn">check</span> '
        '<span class="qc qc-bad">poor</span> &mdash; '
        + "; ".join([rule("pass rate", "pass_rate"), rule("occupancy", "occupancy_pct"),
                     rule("blocking", "blocking_pct")])
        + ". Total row: yield and reads are summed, pass rate is pooled over all reads, "
          "<i>italic</i> values are means across flow cells."
    )


# ---------------------------------------------------------------------------
# INPUT DISCOVERY AND PAIRING
# ---------------------------------------------------------------------------

# MinKNOW file kinds; each file is named <kind>_<run id>.<ext>.
FILE_KINDS = ("pore_activity", "throughput", "sequencing_summary")
FILE_EXTS = {"pore_activity": ".csv", "throughput": ".csv",
             "sequencing_summary": (".txt", ".tsv")}


def _file_kind(path):
    name = Path(path).name
    for kind in FILE_KINDS:
        if name.startswith(kind + "_") and name.endswith(FILE_EXTS[kind]):
            return kind
    return None


def _run_id(path, kind):
    stem = Path(path).stem
    return stem[len(kind) + 1:] if stem.startswith(kind + "_") else stem


def discover_inputs(paths):
    """Expand files and folders (searched recursively) into (kind, path) pairs.

    Returns (items, ignored): items in the order given, folder contents sorted;
    ignored lists files that aren't MinKNOW pore_activity/throughput/summary files.
    """
    items, ignored = [], []
    for p in map(Path, paths):
        if p.is_dir():
            files = sorted(f for f in p.rglob("*") if f.is_file() and _file_kind(f))
            if not files:
                print(f"WARNING: no MinKNOW files found in {p}")
        elif p.is_file():
            files = [p]
        else:
            sys.exit(f"ERROR: {p} does not exist")
        for f in files:
            kind = _file_kind(f)
            if kind:
                items.append((kind, str(f)))
            else:
                ignored.append(str(f))
    return items, ignored


def pair_runs(items):
    """Group (kind, path) items into runs by run ID.

    Returns (runs, unpaired). Runs keep the order in which their first file
    was given; each has run_id, pore_activity, throughput and, if found,
    sequencing_summary. Only renamed files (names not in MinKNOW's
    <kind>_<run id> form, given via --pore-activity etc.) fall back to
    pairing in the order given: a MinKNOW-named file whose run ID matches no
    run belongs to another run and is reported as unpaired, never guessed.
    """
    by_run = {}
    for kind, path in items:
        entry = by_run.setdefault(_run_id(path, kind), {"run_id": _run_id(path, kind)})
        if kind in entry:
            if Path(entry[kind]).resolve() == Path(path).resolve():
                continue  # same file reached twice (e.g. folder and file)
            sys.exit(f"ERROR: two {kind} files for run {entry['run_id']}:\n"
                     f"  {entry[kind]}\n  {path}")
        entry[kind] = path

    runs = [e for e in by_run.values() if "pore_activity" in e and "throughput" in e]
    rest = [e for e in by_run.values() if e not in runs]
    unpaired = [e[k] for e in rest for k in FILE_KINDS if k in e and _file_kind(e[k]) == k]
    spare = {k: [e[k] for e in rest if k in e and _file_kind(e[k]) != k] for k in FILE_KINDS}

    if spare["pore_activity"] and len(spare["pore_activity"]) == len(spare["throughput"]):
        print("NOTE: pairing pore_activity and throughput files without a shared run ID "
              "in the order given.")
        for pa, tp in zip(spare["pore_activity"], spare["throughput"]):
            runs.append({"run_id": _run_id(pa, "pore_activity"), "pore_activity": pa, "throughput": tp})
        spare["pore_activity"], spare["throughput"] = [], []

    # Summaries of complete runs are matched by run ID above; renamed ones only
    # go to runs without one.
    lacking = [r for r in runs if "sequencing_summary" not in r]
    if spare["sequencing_summary"] and len(spare["sequencing_summary"]) == len(lacking):
        print("NOTE: matching sequencing summaries without a shared run ID in the order given.")
        for run, ss in zip(lacking, spare["sequencing_summary"]):
            run["sequencing_summary"] = ss
        spare["sequencing_summary"] = []

    return runs, unpaired + [f for files in spare.values() for f in files]


def assign_labels(runs, labels):
    """Set run["flow_cell"] and run["label"].

    labels: plain labels apply in run order; NAME=LABEL applies to the run
    whose run ID or flow cell ID is NAME. Default: the flow cell ID, numbered
    if a flow cell has several runs.
    """
    mapping = dict(l.split("=", 1) for l in labels if "=" in l)
    ordered = [l for l in labels if "=" not in l]
    if len(ordered) > len(runs):
        print(f"WARNING: {len(ordered)} labels given for {len(runs)} runs; extra labels ignored.")
    flow_cells = [r["run_id"].split("_")[0] for r in runs]
    used = set()
    for i, run in enumerate(runs):
        fc = run["flow_cell"] = flow_cells[i]
        default = fc if flow_cells.count(fc) == 1 else f"{fc} ({flow_cells[:i + 1].count(fc)})"
        label = ordered[i] if i < len(ordered) else default
        for key in (run["run_id"], fc):
            if key in mapping:
                label = mapping[key]
                used.add(key)
                break
        run["label"] = label
    for key in set(mapping) - used:
        print(f"WARNING: label {key}={mapping[key]} matches no run or flow cell ID.")


# ---------------------------------------------------------------------------
# MAIN GENERATION FUNCTION
# ---------------------------------------------------------------------------

def generate_report(runs, title, out_path, sample_hz, nx_values=None):
    """runs: from pair_runs() with labels from assign_labels()."""
    if not runs:
        sys.exit("ERROR: No matching pore_activity/throughput file pairs found.")
    summaries = [r.get("sequencing_summary") for r in runs]
    use_summaries = nx_values is None and all(summaries)
    if nx_values is None and any(summaries) and not all(summaries):
        missing = ", ".join(r["label"] for r in runs if not r.get("sequencing_summary"))
        print(f"WARNING: Nx plot skipped: no sequencing summary for {missing}.")

    print(f"Found {len(runs)} run(s):")

    samples_metrics = []
    tp_data_list = []

    for i, run in enumerate(runs):
        label = run["label"]
        color = SAMPLE_COLORS[i % len(SAMPLE_COLORS)]
        print(f"  [{i+1}] {label}" + (f"  (run {run['run_id']})" if label != run["run_id"] else ""))
        print(f"       pore_activity      : {run['pore_activity']}")
        print(f"       throughput         : {run['throughput']}")
        if run.get("sequencing_summary"):
            print(f"       sequencing_summary : {run['sequencing_summary']}")

        pa_data = parse_pore_activity(run["pore_activity"])
        tp_data = parse_throughput(run["throughput"])
        tp_data_list.append(tp_data)

        metrics = compute_sample_metrics(label, pa_data, tp_data, color, run["flow_cell"])
        samples_metrics.append(metrics)

    # Build JSON payloads
    active_pores_data = build_active_pores_plot_data(samples_metrics, sample_hz)
    yield_data        = build_cumulative_yield_data(samples_metrics, tp_data_list)
    pore_act_data     = build_avg_pore_activity_data(samples_metrics)
    nx_section = ""
    if nx_values is not None:
        try:
            nx_data = build_nx_plot_data_from_values(samples_metrics, nx_values)
        except ValueError as exc:
            sys.exit(f"ERROR: Invalid precomputed Nx values: {exc}")
    elif use_summaries:
        try:
            nx_data = build_nx_plot_data(samples_metrics, summaries)
        except (OSError, ValueError) as exc:
            sys.exit(f"ERROR: Unable to read sequencing summary: {exc}")
    if nx_values is not None or use_summaries:
        nx_iframe = _render_nx_iframe(nx_data, 440)
        nx_section = (
            "<h3>Read Length Nx</h3>\n"
            "<div class=\"plot-container\">\n"
            f"  <iframe srcdoc=\"{nx_iframe}\" height=\"560\" "
            "frameborder=\"0\" scrolling=\"no\"></iframe>\n"
            "</div>\n"
        )

    # Render iframes (escaped for HTML attributes)
    active_pores_iframe  = _render_active_pores_iframe(active_pores_data, 600)
    yield_iframe         = _render_line_chart_iframe(yield_data, 400, "Yield (Gb)")
    pore_activity_iframe = _render_stacked_bar_iframe(pore_act_data)

    summary_rows = _render_summary_rows(samples_metrics)
    total_row    = _render_total_row(samples_metrics)

    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M")
    # d3 goes last: it is large and must not be searched for later tokens.
    report = fill(
        HTML_TEMPLATE,
        title=title,
        table_csv_js=TABLE_CSV_JS,
        summary_header=_render_summary_header(),
        table_note=_render_table_note(),
        summary_rows=summary_rows,
        total_row=total_row,
        active_pores_iframe=active_pores_iframe,
        yield_iframe=yield_iframe,
        nx_section=nx_section,
        pore_activity_iframe=pore_activity_iframe,
        footer=footer_html(title, generated_on),
        page_tip=PAGE_TIP,
        d3_source=D3_SOURCE_TAG,
        d3_js=load_d3(),
    )

    Path(out_path).write_text(report, encoding="utf-8")
    print(f"\nReport written to: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s [INPUT ...] [-o OUT] [-t TITLE] [-l LABEL ...]",
        description="Generate an interactive ONT sequencing HTML report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("inputs", nargs="*", metavar="INPUT",
                        help="MinKNOW files and/or run folders (default: current folder)")
    parser.add_argument("-o", "--out", default="ont-report.html",
                        help="Output HTML file (default: ont-report.html)")
    parser.add_argument("-t", "--title", default="ONT Sequencing Report",
                        help="Report title (default: 'ONT Sequencing Report')")
    parser.add_argument("-l", "--labels", nargs="+", metavar="LABEL", default=[],
                        help="NAME=LABEL (NAME = flow cell or run ID), or labels in run order")
    parser.add_argument("--sample-hz", type=int, default=5, metavar="MIN",
                        help="Sampling interval for the active-pores plot in minutes (default: 5)")
    explicit = parser.add_argument_group("files of one kind (e.g. renamed files)")
    explicit.add_argument("--pore-activity", nargs="+", metavar="FILE", default=[])
    explicit.add_argument("--throughput", nargs="+", metavar="FILE", default=[])
    explicit.add_argument("--sequencing-summary", nargs="+", metavar="FILE", default=[])
    parser.add_argument("--nx-values", metavar="JSON",
                        help="Precomputed Nx values as JSON, one array of {x,y} points per run")
    args = parser.parse_args()
    try:
        nx_values = json.loads(args.nx_values) if args.nx_values else None
    except json.JSONDecodeError as exc:
        parser.error(f"--nx-values must be valid JSON: {exc}")

    has_explicit = args.pore_activity or args.throughput or args.sequencing_summary
    items, ignored = discover_inputs(args.inputs or ([] if has_explicit else ["."]))
    for kind in FILE_KINDS:
        items += [(kind, f) for f in getattr(args, kind)]
    if ignored:
        print("NOTE: ignoring files that aren't MinKNOW pore_activity/throughput/"
              "sequencing_summary files:\n  " + "\n  ".join(ignored))

    runs, unpaired = pair_runs(items)
    if unpaired:
        print("WARNING: these files match no complete run (pore_activity + throughput) "
              "and are ignored:\n  " + "\n  ".join(unpaired))
    if not runs:
        sys.exit("ERROR: No pore_activity/throughput pairs found. Pass MinKNOW run folders "
                 "or files, e.g.: generate-ont-report.py /data/run1 /data/run2")
    assign_labels(runs, args.labels)

    generate_report(runs, title=args.title, out_path=args.out,
                    sample_hz=args.sample_hz, nx_values=nx_values)


if __name__ == "__main__":
    main()
