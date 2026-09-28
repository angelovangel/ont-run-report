#!/usr/bin/env python3
"""
generate-hifi-report.py
========================
Generate an interactive HTML report for PacBio HiFi (Revio/Sequel) runs
from one or more "run-qc-export.csv" files, where each row of a CSV
describes a single SMRT Cell.

Usage
-----
  python generate-hifi-report.py FILE [FILE ...] [OPTIONS]

Options
-------
  --title   TEXT   Report title (default: "HiFi Sequencing Report").
  --out     PATH   Output HTML file (default: report.html).

If no FILE arguments are given, the script auto-discovers files matching
"*qc-export*.csv" in the current directory.

Examples
--------
  # Explicit files
  python generate-hifi-report.py 20260913_P12U8-run-qc-export.csv \\
                                  20260914_P12U8-run-qc-export.csv \\
                                  --title "P12U8 HiFi Report" --out p12u8.html

  # Auto-discover all *qc-export*.csv in the current directory
  python generate-hifi-report.py
"""

import argparse
import csv
import glob
import json
import re
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# COLOUR PALETTE (one colour per unique Sample Name)
# ---------------------------------------------------------------------------
SAMPLE_COLORS = [
    "#E74C3C", "#3498DB", "#2ECC71", "#F39C12",
    "#9B59B6", "#1ABC9C", "#E67E22", "#34495E",
    "#E91E63", "#00BCD4",
]


# ---------------------------------------------------------------------------
# VALUE PARSING HELPERS
# ---------------------------------------------------------------------------

def _clean(v):
    if v is None:
        return None
    v = str(v).strip()
    return v if v else None


def _float(v):
    v = _clean(v)
    if v is None:
        return None
    try:
        return float(v.replace(",", ""))
    except ValueError:
        return None


def _int(v):
    f = _float(v)
    return int(round(f)) if f is not None else None


def _frac_to_pct(v):
    """Some exported '(%)' columns are actually stored as a 0-1 fraction."""
    f = _float(v)
    return round(f * 100, 2) if f is not None else None


def _q_value(v):
    """'Read quality (median)' is stored like 'Q34' -> numeric 34."""
    v = _clean(v)
    if v is None:
        return None
    m = re.match(r"Q?\s*(\d+(?:\.\d+)?)", v)
    return float(m.group(1)) if m else None


def _parse_dt(v):
    v = _clean(v)
    if v is None:
        return None
    for fmt in ("%m.%d.%Y %H:%M", "%m/%d/%Y %H:%M", "%Y-%m-%d %H:%M"):
        try:
            return datetime.strptime(v, fmt)
        except ValueError:
            continue
    return None


# ---------------------------------------------------------------------------
# METRIC DEFINITIONS
#   key      -- identifier used internally / in JS
#   csv_col  -- exact column name in the run-qc-export.csv
#   label    -- human-readable label (table header / dropdown option)
#   parse    -- parsing function
#   decimals -- rounding for float display (None for ints)
#   in_table -- include this metric as a column in the main summary table
#               (all metrics are always available in the plot dropdowns)
# ---------------------------------------------------------------------------
MASTER_METRICS = [
    ("hifi_yield_gb",             "HiFi Yield (Gb)",              "HiFi Yield (Gb)",            _float,      2, True),
    ("total_bases_gb",            "Total Bases (Gb)",             "Total Bases (Gb)",           _float,      2, False),
    ("polymerase_rl_bp",          "Polymerase RL (bp)",           "Polymerase RL (bp)",         _int,        0, True),
    ("polymerase_n50_bp",         "Polymerase N50 (bp)",          "Polymerase N50 (bp)",        _int,        0, True),
    ("longest_subread_bp",        "Longest Subread (bp)",         "Longest Subread (bp)",       _int,        0, False),
    ("longest_subread_n50_bp",    "Longest Subread N50 (bp)",     "Longest Subread N50 (bp)",   _int,        0, False),
    ("q20_reads",                 "Q20+ reads",                   "Q20+ Reads",                 _int,        0, False),
    ("mean_length_bp",            "Mean Length",                  "Mean Length (bp)",           _int,        0, True),
    ("read_quality_q",            "Read quality (median)",        "Read Quality (median, Q)",   _q_value,    0, True),
    ("q30_bq_pct",                "Q30+ BQ (%)",                  "Q30+ BQ (%)",                _frac_to_pct, 2, True),
    ("loading_conc_pM",           "Loading concentration (pM)",   "Loading Conc. (pM)",         _float,      0, False),
    ("p0_pct",                    "P0 %",                         "P0 (%)",                     _float,      2, False),
    ("p1_pct",                    "P1 %",                         "P1 (%)",                     _float,      2, True),
    ("p2_pct",                    "P2 %",                         "P2 (%)",                     _float,      2, False),
    ("local_base_rate",           "Local Base Rate",              "Local Base Rate (bp/s)",     _float,      2, False),
    ("missing_adapter_pct",       "Missing Adapter (%)",          "Missing Adapter (%)",        _float,      2, False),
    ("control_concordance_mean",  "Control Concordance Mean",     "Control Concordance Mean",   _float,      4, False),
    ("insert_size_bp",            "Insert Size (bp)",             "Insert Size (bp)",           _int,        0, False),
    ("movie_time_hours",          "Movie Time (hours)",           "Movie Time (h)",             _float,      1, False),
    ("acquisition_no",            "Use",                          "Acquisition #",              _int,        0, False),
]

# Table columns shown in the main summary table (subset of MASTER_METRICS,
# by key, in display order). All metrics remain available in the plot
# dropdowns regardless of this list.
TABLE_METRIC_KEYS = [
    "hifi_yield_gb", "q20_reads", "polymerase_rl_bp", "polymerase_n50_bp",
    "mean_length_bp", "read_quality_q", "q30_bq_pct", "p1_pct",
]


# ---------------------------------------------------------------------------
# CSV PARSING
# ---------------------------------------------------------------------------

def parse_qc_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        return [row for row in reader]


def process_row(raw, source_file):
    rec = {
        "run_name":      (raw.get("Run Name") or "").strip(),
        "run_id":        (raw.get("Run ID") or "").strip(),
        "sample_name":   (raw.get("Sample Name") or "").strip(),
        "well":          (raw.get("Sample Well") or "").strip(),
        "cell_id":       (raw.get("Cell ID") or "").strip(),
        "status":        (raw.get("Status") or "").strip(),
        "instrument_sn": (raw.get("Instrument SN") or "").strip(),
        "run_start_raw": (raw.get("Run Start") or "").strip(),
        "source_file":   source_file,
    }
    rec["run_start_dt"] = _parse_dt(rec["run_start_raw"])
    rec["acquisition_no"] = _int(raw.get("Use"))
    rec["multi_use"] = (raw.get("Is Multi-Use Cell") or "").strip().upper() == "TRUE"

    metrics = {}
    for key, col, _label, fn, dec, _table in MASTER_METRICS:
        val = fn(raw.get(col, ""))
        if val is not None and dec is not None and isinstance(val, float):
            val = round(val, dec)
        metrics[key] = val
    rec["metrics"] = metrics

    date_part = rec["run_name"].split("_")[0] if rec["run_name"] else ""
    rec["label"] = f"{date_part} {rec['well']}".strip() or rec["cell_id"] or "cell"
    return rec


def load_all_cells(paths):
    cells = []
    seen = set()
    for path in paths:
        rows = parse_qc_csv(path)
        for raw in rows:
            rec = process_row(raw, Path(path).name)
            dedup_key = (rec["run_id"], rec["well"]) if rec["run_id"] else None
            if dedup_key and dedup_key in seen:
                continue
            if dedup_key:
                seen.add(dedup_key)
            cells.append(rec)

    cells.sort(key=lambda r: (r["run_start_dt"] or datetime.min, r["well"]))
    return cells


# Color-by grouping modes, in dropdown display order. Key is used internally
# / in JS; label is the human-readable dropdown option text.
COLOR_GROUP_MODES = [
    ("sample", "Sample"),
    ("run", "Run"),
    ("cell", "Flow Cell"),
]


def _group_value(rec, mode):
    """Value used to group/color a cell for a given color-by mode."""
    if mode == "sample":
        return rec["sample_name"] or rec["run_name"] or "Unknown"
    if mode == "run":
        return rec["run_name"] or "Unknown"
    if mode == "cell":
        return rec["cell_id"] or "Unknown"
    return "Unknown"


def assign_colors(cells):
    """Assign each cell a color for every color-by mode (sample/run/cell).

    Populates rec["colors"] = {"sample": "#..", "run": "#..", "cell": "#.."}
    and rec["color"] as a backward-compatible alias for the sample mode.
    Returns {mode: {group_name: color}} for each mode.
    """
    color_maps = {mode: {} for mode, _label in COLOR_GROUP_MODES}
    for rec in cells:
        colors = {}
        for mode, _label in COLOR_GROUP_MODES:
            name = _group_value(rec, mode)
            cmap = color_maps[mode]
            if name not in cmap:
                cmap[name] = SAMPLE_COLORS[len(cmap) % len(SAMPLE_COLORS)]
            colors[mode] = cmap[name]
        rec["colors"] = colors
        rec["color"] = colors["sample"]
    return color_maps


# ---------------------------------------------------------------------------
# TABLE RENDERING
# ---------------------------------------------------------------------------

METRICS_BY_KEY = {key: (label, dec) for key, _col, label, _fn, dec, _t in MASTER_METRICS}

# Unified column model: identity columns (text, some filterable via a
# dropdown), metric columns (numeric, pulled from rec["metrics"]), and the
# derived "Acquisition" column (Use + Is Multi-Use Cell combined).
COLUMNS = [
    {"key": "run_name",    "label": "Run Name", "type": "text", "filter": True},
    {"key": "cell_id",     "label": "Cell ID",   "type": "text", "filter": True},
    {"key": "acquisition", "label": "Acquisition", "type": "number",
     "dec": 0, "metric": False, "filter": False},
    {"key": "sample_name", "label": "Sample",    "type": "text", "filter": True},
    {"key": "well",        "label": "Well",      "type": "text", "filter": False},
]
for _mkey in TABLE_METRIC_KEYS:
    _label, _dec = METRICS_BY_KEY[_mkey]
    COLUMNS.append({"key": _mkey, "label": _label, "type": "number",
                     "dec": _dec, "metric": True, "filter": False})

FILTER_COLUMNS = [c for c in COLUMNS if c["filter"]]


def _fmt_number(dec, value):
    if value is None:
        return "-"
    return f"{value:,.0f}" if dec == 0 else f"{value:,.{dec}f}"


def _html_escape(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
                  .replace(">", "&gt;").replace('"', "&quot;"))


def _col_value(rec, col):
    key = col["key"]
    if key in ("run_name", "cell_id", "sample_name", "well"):
        return rec[key]
    if key == "acquisition":
        return rec["acquisition_no"]
    return rec["metrics"].get(key)


def _col_display(rec, col):
    val = _col_value(rec, col)
    if col["type"] == "text":
        return val or ""
    if col["key"] == "acquisition":
        if val is None:
            return "-"
        return f"{val:,.0f} (multi-use)" if rec["multi_use"] else f"{val:,.0f}"
    return _fmt_number(col["dec"], val)


def render_table_header():
    ths = []
    for col in COLUMNS:
        align = ' style="text-align:right"' if col["type"] == "number" else ""
        is_metric = " data-metric=\"true\"" if col.get("metric") else ""
        ths.append(f'<th{align} data-type="{col["type"]}"{is_metric} onclick="sortTable(this)">'
                   f'{col["label"]}<span class="sort-arrow"></span></th>')
    return "".join(ths)


def render_table_row(rec):
    tds = []
    for i, col in enumerate(COLUMNS):
        val = _col_value(rec, col)
        disp = _col_display(rec, col)
        if col["type"] == "number":
            data_value = "" if val is None else repr(val)
        else:
            data_value = (val or "").lower()

        style_parts = []
        class_parts = []
        if col["type"] == "number":
            style_parts.append("text-align:right")
        if i == 0:
            class_parts.append("color-stripe")
            style_parts.append(f"border-left:4px solid {rec['colors']['sample']}")
        style_attr = f' style="{";".join(style_parts)}"' if style_parts else ""
        class_attr = f' class="{" ".join(class_parts)}"' if class_parts else ""

        tds.append(f'<td{class_attr}{style_attr} data-value="{_html_escape(data_value)}">'
                   f'{_html_escape(disp)}</td>')

    row_attrs = (f'data-run="{_html_escape(rec["run_name"])}" '
                f'data-cell="{_html_escape(rec["cell_id"])}" '
                f'data-sample="{_html_escape(rec["sample_name"])}" '
                f'data-color-sample="{rec["colors"]["sample"]}" '
                f'data-color-run="{rec["colors"]["run"]}" '
                f'data-color-cell="{rec["colors"]["cell"]}"')
    return f"    <tr {row_attrs}>" + "".join(tds) + "</tr>"


def render_filter_options(cells, key):
    values = sorted({rec[key] for rec in cells if rec[key]})
    return "".join(f'<option value="{_html_escape(v)}">{_html_escape(v)}</option>'
                   for v in values)


def render_table_footer(cells):
    """Build a <tfoot> total row summing HiFi Yield and HiFi Reads."""
    n_text_cols = sum(1 for c in COLUMNS if c["type"] == "text" or c["key"] == "acquisition")
    total_yield = sum(
        (c["metrics"].get("hifi_yield_gb") or 0) for c in cells
    )
    total_reads = sum(
        int(c["metrics"].get("q20_reads") or 0) for c in cells
    )
    tds = []
    for i, col in enumerate(COLUMNS):
        style = ' style="text-align:right"' if col["type"] == "number" else ""
        if i == 0:
            tds.append(f"<td><strong>Total</strong></td>")
        elif col["key"] == "hifi_yield_gb":
            tds.append(f"<td style='text-align:right'>{total_yield:,.2f}</td>")
        elif col["key"] == "q20_reads":
            tds.append(f"<td style='text-align:right'>{total_reads:,}</td>")
        else:
            tds.append(f"<td{style}></td>")
    return f"    <tr>" + "".join(tds) + "</tr>"


# ---------------------------------------------------------------------------
# CHART DATA
# ---------------------------------------------------------------------------

def build_chart_payload(cells):
    metrics_meta = []
    for key, _col, label, _fn, dec, _t in MASTER_METRICS:
        if any(rec["metrics"].get(key) is not None for rec in cells):
            metrics_meta.append({"key": key, "label": label, "decimals": dec})

    records = []
    for rec in cells:
        records.append({
            "label": rec["label"],
            "sample": rec["sample_name"] or rec["run_name"],
            "well": rec["well"],
            "run": rec["run_name"],
            "cellId": rec["cell_id"],
            "status": rec["status"],
            "colors": rec["colors"],
            "metrics": rec["metrics"],
        })

    # Legend entries per color-by mode: [{name, color}, ...]
    groups = {}
    for mode, _label in COLOR_GROUP_MODES:
        seen = {}
        for rec in cells:
            name = _group_value(rec, mode)
            if name not in seen:
                seen[name] = rec["colors"][mode]
        groups[mode] = [{"name": n, "color": c} for n, c in seen.items()]

    color_modes = [{"key": key, "label": label} for key, label in COLOR_GROUP_MODES]

    return {"records": records, "metrics": metrics_meta,
            "groups": groups, "colorModes": color_modes}


# ---------------------------------------------------------------------------
# HTML TEMPLATE (outer page)
# ---------------------------------------------------------------------------

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<meta name="description" content="Interactive PacBio HiFi run QC report.">
<style>
  *,*::before,*::after{{box-sizing:border-box;}}
  body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
        max-width:1100px;margin:0 auto;padding:24px 28px;color:#2d3748;
        line-height:1.55;background:#f7f8fa;}}
  h2{{color:#1a202c;border-bottom:2px solid #e2e8f0;padding-bottom:8px;margin-top:40px;}}
  h3{{color:#2d3748;margin-top:28px;font-size:1.05rem;}}
  .subtitle{{color:#718096;font-size:13.5px;margin-top:-8px;}}
  .table-wrap{{overflow-x:auto;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.07);}}
  table{{border-collapse:collapse;width:100%;margin:14px 0;font-size:13px;
         background:#fff;white-space:nowrap;}}
  th{{background:#edf2f7;color:#2d3748;font-weight:600;text-align:left;
      padding:9px 11px;border-bottom:2px solid #cbd5e0;position:sticky;top:0;
      cursor:pointer;user-select:none;}}
  th:hover{{background:#e2e8f0;}}
  .sort-arrow{{font-size:10px;color:#4a5568;margin-left:3px;}}
  td{{padding:8px 13px 14px;border-bottom:1px solid #e9ecef;position:relative;vertical-align:middle;}}
  tr:last-child td{{border-bottom:none;}}
  tr:hover td{{background:rgba(59,130,246,.04);}}
  .bar-track{{position:absolute;bottom:4px;left:10px;right:10px;height:3px;
              background:rgba(59,130,246,.12);border-radius:99px;overflow:hidden;}}
  .bar-fill{{position:absolute;inset:0;height:100%;border-radius:99px;
             background:#3b82f6;
             transition:width .35s cubic-bezier(.4,0,.2,1);}}
  tfoot tr td{{background:#edf2f7;font-weight:600;border-top:2px solid #cbd5e0;
               color:#2d3748;padding:8px 11px;white-space:nowrap;}}
  .filter-bar{{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:14px 0 4px;}}
  .filter-bar label{{font-size:12.5px;font-weight:600;color:#555;}}
  .filter-bar select{{padding:5px 10px;font-size:12.5px;border:1px solid #cbd5e0;
                      border-radius:6px;background:#f8fafc;color:#2d3748;cursor:pointer;}}
  .row-count{{font-size:12.5px;color:#718096;margin-left:auto;}}
  .plot-container{{margin:20px 0;background:#fff;border-radius:8px;
                   box-shadow:0 1px 4px rgba(0,0,0,.08);border:1px solid #e2e8f0;overflow:hidden;}}
  .footer{{margin-top:48px;border-top:1px solid #e2e8f0;padding-top:14px;
           color:#a0aec0;font-size:12px;text-align:center;}}
</style>
</head>
<body>

<h2>{title}</h2>
<p class="subtitle">{n_cells} SMRT cell(s) across {n_runs} run(s) &mdash; generated on {generated_on}</p>

<h3>Cell Summary</h3>
<div class="filter-bar">
  <label for="filter-run">Run Name</label>
  <select id="filter-run"><option value="">All</option>{filter_run_options}</select>
  <label for="filter-cell">Cell ID</label>
  <select id="filter-cell"><option value="">All</option>{filter_cell_options}</select>
  <label for="filter-sample">Sample</label>
  <select id="filter-sample"><option value="">All</option>{filter_sample_options}</select>
  <label for="color-by-select">Color by</label>
  <select id="color-by-select">
    <option value="sample">Sample</option>
    <option value="run">Run</option>
    <option value="cell">Flow Cell</option>
  </select>
  <span class="row-count" id="row-count"></span>
</div>
<div id="table-legend" style="display:flex;flex-wrap:wrap;gap:12px;margin:2px 0 12px;font-size:12px;color:#555;"></div>
<div class="table-wrap">
<table id="cell-table">
  <thead><tr>{table_header}</tr></thead>
  <tbody id="cell-tbody">
{table_rows}
  </tbody>
  <tfoot>
{table_footer}
  </tfoot>
</table>
</div>
<script>
  function sortTable(th) {{
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    const headerCells = Array.from(th.parentNode.children);
    const colIndex = headerCells.indexOf(th);
    const type = th.dataset.type;
    const dir = th.dataset.dir === 'asc' ? 'desc' : 'asc';
    headerCells.forEach(h => {{ h.dataset.dir = ''; h.querySelector('.sort-arrow').textContent = ''; }});
    th.dataset.dir = dir;
    th.querySelector('.sort-arrow').textContent = dir === 'asc' ? '\u25B2' : '\u25BC';

    const rows = Array.from(tbody.querySelectorAll('tr'));
    rows.sort((a, b) => {{
      const av = a.children[colIndex].dataset.value;
      const bv = b.children[colIndex].dataset.value;
      let cmp;
      if (type === 'number') {{
        const an = av === '' ? -Infinity : parseFloat(av);
        const bn = bv === '' ? -Infinity : parseFloat(bv);
        cmp = an - bn;
      }} else {{
        cmp = av.localeCompare(bv);
      }}
      return dir === 'asc' ? cmp : -cmp;
    }});
    rows.forEach(r => tbody.appendChild(r));
  }}

  function applyFilters() {{
    const run = document.getElementById('filter-run').value;
    const cell = document.getElementById('filter-cell').value;
    const sample = document.getElementById('filter-sample').value;
    const rows = document.querySelectorAll('#cell-tbody tr');
    let shown = 0;
    rows.forEach(r => {{
      const match = (!run || r.dataset.run === run) &&
                    (!cell || r.dataset.cell === cell) &&
                    (!sample || r.dataset.sample === sample);
      r.style.display = match ? '' : 'none';
      if (match) shown++;
    }});
    document.getElementById('row-count').textContent = `Showing ${{shown}} of ${{rows.length}} cell(s)`;
  }}

  ['filter-run', 'filter-cell', 'filter-sample'].forEach(id =>
    document.getElementById(id).addEventListener('change', applyFilters));
  applyFilters();

  const COLOR_ATTR = {{ sample: 'colorSample', run: 'colorRun', cell: 'colorCell' }};
  const NAME_ATTR  = {{ sample: 'sample',       run: 'run',       cell: 'cell' }};

  function applyRowColors(mode) {{
    document.querySelectorAll('#cell-tbody tr').forEach(r => {{
      const stripe = r.querySelector('td.color-stripe');
      if (stripe) stripe.style.borderLeftColor = r.dataset[COLOR_ATTR[mode]];
    }});
  }}

  function buildTableLegend(mode) {{
    const seen = new Map();
    document.querySelectorAll('#cell-tbody tr').forEach(r => {{
      const name = r.dataset[NAME_ATTR[mode]] || '(none)';
      if (!seen.has(name)) seen.set(name, r.dataset[COLOR_ATTR[mode]]);
    }});
    const legend = document.getElementById('table-legend');
    legend.innerHTML = '';
    seen.forEach((color, name) => {{
      const item = document.createElement('div');
      item.style.display = 'flex';
      item.style.alignItems = 'center';
      item.style.gap = '5px';
      item.innerHTML = `<span style="width:10px;height:10px;border-radius:3px;` +
        `background:${{color}};display:inline-block;"></span>${{name}}`;
      legend.appendChild(item);
    }});
  }}

  document.getElementById('color-by-select').addEventListener('change', e => {{
    applyRowColors(e.target.value);
    buildTableLegend(e.target.value);
  }});
  applyRowColors('sample');
  buildTableLegend('sample');

  // ---------------------------------------------------------------------------
  // DATA BARS: thin inset pill at the bottom of each numeric metric cell
  // Largest visible value per column = 100% fill. Acquisition excluded.
  // ---------------------------------------------------------------------------
  const METRIC_COL_INDICES = (function() {{
    const ths = Array.from(document.querySelectorAll('#cell-table thead th'));
    return ths.map((th, i) => th.dataset.metric === 'true' ? i : -1).filter(i => i >= 0);
  }})();

  function ensureTrack(td) {{
    let track = td.querySelector('.bar-track');
    if (!track) {{
      track = document.createElement('div');
      track.className = 'bar-track';
      const fill = document.createElement('div');
      fill.className = 'bar-fill';
      track.appendChild(fill);
      td.appendChild(track);
    }}
    return track.querySelector('.bar-fill');
  }}

  function applyDataBars() {{
    const rows = Array.from(
      document.querySelectorAll('#cell-tbody tr')
    ).filter(r => r.style.display !== 'none');

    METRIC_COL_INDICES.forEach(ci => {{
      let maxVal = 0;
      rows.forEach(r => {{
        const v = parseFloat(r.children[ci]?.dataset.value);
        if (!isNaN(v) && v > maxVal) maxVal = v;
      }});
      rows.forEach(r => {{
        const td = r.children[ci];
        if (!td) return;
        const fill = ensureTrack(td);
        const v = parseFloat(td.dataset.value);
        fill.style.width = (!isNaN(v) && maxVal > 0)
          ? `${{Math.round((v / maxVal) * 100)}}%`
          : '0%';
      }});
    }});
  }}

  // Patch sortTable to re-apply bars after each sort
  (function() {{
    const _orig = sortTable;
    window.sortTable = function(th) {{ _orig(th); applyDataBars(); }};
  }})();

  // Re-apply bars when filters change
  ['filter-run', 'filter-cell', 'filter-sample'].forEach(id =>
    document.getElementById(id).addEventListener('change', applyDataBars));

  applyDataBars();
</script>

<h3>Metric by Cell</h3>
<div class="plot-container">
  <iframe srcdoc="{bar_iframe}"
          width="100%" height="600" frameborder="0"
          scrolling="no" style="border:none;display:block;"></iframe>
</div>

<h3>Metric Correlation</h3>
<div class="plot-container">
  <iframe srcdoc="{scatter_iframe}"
          width="100%" height="600" frameborder="0"
          scrolling="no" style="border:none;display:block;"></iframe>
</div>

<div class="footer">
  Generated on {generated_on} &mdash; {title}
</div>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# IFRAME: shared CSS block
# ---------------------------------------------------------------------------

_IFRAME_CSS = """\
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:#fff;margin:0;padding:20px;color:#333;}}
    #chart-container{{max-width:1020px;margin:0 auto;position:relative;}}
    .tooltip{{position:absolute;padding:10px;font-size:13px;
              background:rgba(255,255,255,.95);border:1px solid #ddd;border-radius:4px;
              pointer-events:none;opacity:0;box-shadow:0 4px 6px rgba(0,0,0,.1);
              white-space:nowrap;transition:opacity .2s;z-index:10;}}
    .grid line{{stroke:#ddd;stroke-dasharray:4,4;}}
    .grid path{{stroke-width:0;}}
    .axis text{{font-size:12px;fill:#555;}}
    .axis path,.axis line{{stroke:#bbb;}}
    h2{{text-align:center;color:#222;margin-bottom:6px;font-size:1.1rem;}}
    .controls{{display:flex;justify-content:center;align-items:center;gap:10px;
               margin-bottom:14px;flex-wrap:wrap;}}
    .controls label{{font-size:12.5px;color:#555;font-weight:600;}}
    .controls select{{
      padding:5px 10px;font-size:12.5px;border:1px solid #cbd5e0;border-radius:6px;
      background:#f8fafc;color:#2d3748;cursor:pointer;
    }}
    .legend-text{{font-size:12px;user-select:none;}}
    .bar{{opacity:.9;}}
    .bar:hover{{opacity:1;}}
    .dot{{opacity:.85;stroke:#fff;stroke-width:1px;}}
    .dot:hover{{opacity:1;}}
"""


# ---------------------------------------------------------------------------
# IFRAME 1: BAR CHART -- dropdown selects which metric to plot per cell
# ---------------------------------------------------------------------------

BAR_CHART_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Metric by Cell</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
""" + _IFRAME_CSS + """\
  </style>
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
    <div class="controls">
      <label for="metric-select">Metric</label>
      <select id="metric-select"></select>
      <label for="color-by-select">Color by</label>
      <select id="color-by-select"></select>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const payload = {plot_data_json};
    const records = payload.records;
    const metrics = payload.metrics;
    const groups = payload.groups;
    const colorModes = payload.colorModes;
    let colorMode = colorModes.length ? colorModes[0].key : "sample";

    const select = d3.select("#metric-select");
    select.selectAll("option")
      .data(metrics)
      .join("option")
      .attr("value", d => d.key)
      .text(d => d.label);

    const colorSelect = d3.select("#color-by-select");
    colorSelect.selectAll("option")
      .data(colorModes)
      .join("option")
      .attr("value", d => d.key)
      .text(d => d.label);
    colorSelect.property("value", colorMode);

    const margin = {{top:20,right:160,bottom:70,left:70}};
    const W = 960 - margin.left - margin.right;
    const H = 460 - margin.top - margin.bottom;

    const svg = d3.select("#chart").append("svg")
        .attr("width",  W + margin.left + margin.right)
        .attr("height", H + margin.top  + margin.bottom)
      .append("g")
        .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const tooltip = d3.select("#tooltip");

    const x = d3.scaleBand().range([0, W]).padding(0.3);
    const y = d3.scaleLinear().range([H, 0]);

    const gGrid  = svg.append("g").attr("class","grid");
    const gXAxis = svg.append("g").attr("class","axis").attr("transform",`translate(0,${{H}})`);
    const gYAxis = svg.append("g").attr("class","axis");

    const yLabel = svg.append("text")
       .attr("transform","rotate(-90)")
       .attr("x", -H/2).attr("y", -55)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555");

    const barsG = svg.append("g");

    // Legend -- top-right of the plot area, inside the chart. Click to filter.
    const legendG = svg.append("g").attr("transform", `translate(${{W+14}},0)`);
    let hiddenGroups = new Set();

    function groupNameOf(r) {{
      if (colorMode === "run") return r.run || "Unknown";
      if (colorMode === "cell") return r.cellId || "Unknown";
      return r.sample || "Unknown";
    }}

    function renderLegend() {{
      legendG.selectAll("*").remove();
      (groups[colorMode] || []).forEach((g, i) => {{
        const hidden = hiddenGroups.has(g.name);
        const row = legendG.append("g")
            .attr("transform", `translate(0,${{i*20}})`)
            .style("cursor","pointer")
            .on("click", function() {{
              if (hiddenGroups.has(g.name)) hiddenGroups.delete(g.name);
              else hiddenGroups.add(g.name);
              renderLegend();
              draw(select.property("value"));
            }});
        row.append("rect")
            .attr("width",11).attr("height",11).attr("rx",2)
            .attr("fill", hidden ? "#ccc" : g.color);
        row.append("text")
            .attr("class","legend-text")
            .attr("x",16).attr("y",9)
            .attr("fill", hidden ? "#999" : "#333")
            .style("text-decoration", hidden ? "line-through" : "none")
            .text(g.name);
      }});
    }}

    function draw(metricKey) {{
      const meta = metrics.find(m => m.key === metricKey) || {{}};
      document.getElementById("chart-title").textContent = meta.label || metricKey;

      const data = records.filter(r => r.metrics[metricKey] !== null &&
                                        r.metrics[metricKey] !== undefined &&
                                        !hiddenGroups.has(groupNameOf(r)));

      x.domain(data.map(d => d.label));
      const yMax = d3.max(data, d => d.metrics[metricKey]) || 1;
      y.domain([0, yMax * 1.08]);

      gGrid.call(d3.axisLeft(y).tickSize(-W).tickFormat(""));
      gXAxis.call(d3.axisBottom(x))
            .selectAll("text")
            .attr("transform", "rotate(-30)")
            .style("text-anchor", "end");
      gYAxis.call(d3.axisLeft(y).tickFormat(d3.format(",.2~f")));
      yLabel.text(meta.label || metricKey);

      barsG.selectAll("rect").data(data, d => d.label)
        .join("rect")
        .attr("class", "bar")
        .attr("x", d => x(d.label))
        .attr("width", x.bandwidth())
        .attr("fill", d => d.colors[colorMode])
        .attr("y", d => y(d.metrics[metricKey]))
        .attr("height", d => H - y(d.metrics[metricKey]))
        .on("mousemove", function(event, d) {{
          const val = d.metrics[metricKey];
          const decimals = meta.decimals ?? 2;
          tooltip.transition().duration(40).style("opacity",1);
          tooltip.html(
            `<strong>${{d.label}}</strong><br/>` +
            `Sample: ${{d.sample}}<br/>` +
            `Run: ${{d.run}}<br/>` +
            `Cell ID: ${{d.cellId}}<br/>` +
            `<span style="color:${{d.colors[colorMode]}}">&#9632;</span> ${{meta.label}}: ` +
            `<b>${{d3.format(",." + decimals + "f")(val)}}</b>`
          ).style("left",(event.pageX+15)+"px")
           .style("top",Math.min(event.pageY-28, window.innerHeight-160)+"px");
        }})
        .on("mouseout", () => tooltip.transition().duration(200).style("opacity",0));
    }}

    select.on("change", function() {{ draw(this.value); }});
    colorSelect.on("change", function() {{
      colorMode = this.value;
      hiddenGroups = new Set();
      renderLegend();
      draw(select.property("value"));
    }});
    if (metrics.length) {{
      select.property("value", metrics[0].key);
      renderLegend();
      draw(metrics[0].key);
    }}
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# IFRAME 2: SCATTER -- two dropdowns select X and Y metrics to correlate
# ---------------------------------------------------------------------------

SCATTER_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Metric Correlation</title>
  <script src="https://d3js.org/d3.v7.min.js"></script>
  <style>
""" + _IFRAME_CSS + """\
  </style>
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
    <div class="controls">
      <label for="x-select">X axis</label>
      <select id="x-select"></select>
      <label for="y-select">Y axis</label>
      <select id="y-select"></select>
      <label for="color-by-select">Color by</label>
      <select id="color-by-select"></select>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const payload = {plot_data_json};
    const records = payload.records;
    const metrics = payload.metrics;
    const groups = payload.groups;
    const colorModes = payload.colorModes;
    let colorMode = colorModes.length ? colorModes[0].key : "sample";

    const xSelect = d3.select("#x-select");
    const ySelect = d3.select("#y-select");
    [xSelect, ySelect].forEach(sel => {{
      sel.selectAll("option")
        .data(metrics)
        .join("option")
        .attr("value", d => d.key)
        .text(d => d.label);
    }});

    const colorSelect = d3.select("#color-by-select");
    colorSelect.selectAll("option")
      .data(colorModes)
      .join("option")
      .attr("value", d => d.key)
      .text(d => d.label);
    colorSelect.property("value", colorMode);

    const margin = {{top:20,right:160,bottom:60,left:70}};
    const W = 960 - margin.left - margin.right;
    const H = 460 - margin.top - margin.bottom;

    const svg = d3.select("#chart").append("svg")
        .attr("width",  W + margin.left + margin.right)
        .attr("height", H + margin.top  + margin.bottom)
      .append("g")
        .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const tooltip = d3.select("#tooltip");

    const x = d3.scaleLinear().range([0, W]);
    const y = d3.scaleLinear().range([H, 0]);

    const gGrid  = svg.append("g").attr("class","grid");
    const gXAxis = svg.append("g").attr("class","axis").attr("transform",`translate(0,${{H}})`);
    const gYAxis = svg.append("g").attr("class","axis");

    const xLabel = svg.append("text")
       .attr("x", W/2).attr("y", H+45)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555");
    const yLabel = svg.append("text")
       .attr("transform","rotate(-90)")
       .attr("x", -H/2).attr("y", -55)
       .attr("text-anchor","middle").attr("font-size","12px").attr("fill","#555");

    const dotsG = svg.append("g");

    // Legend -- top-right of the plot area, inside the chart. Click to filter.
    const legendG = svg.append("g").attr("transform", `translate(${{W+14}},0)`);
    let hiddenGroups = new Set();

    function groupNameOf(r) {{
      if (colorMode === "run") return r.run || "Unknown";
      if (colorMode === "cell") return r.cellId || "Unknown";
      return r.sample || "Unknown";
    }}

    function renderLegend() {{
      legendG.selectAll("*").remove();
      (groups[colorMode] || []).forEach((g, i) => {{
        const hidden = hiddenGroups.has(g.name);
        const row = legendG.append("g")
            .attr("transform", `translate(0,${{i*20}})`)
            .style("cursor","pointer")
            .on("click", function() {{
              if (hiddenGroups.has(g.name)) hiddenGroups.delete(g.name);
              else hiddenGroups.add(g.name);
              renderLegend();
              draw();
            }});
        row.append("circle")
            .attr("cx",5.5).attr("cy",5.5).attr("r",5.5)
            .attr("fill", hidden ? "#ccc" : g.color);
        row.append("text")
            .attr("class","legend-text")
            .attr("x",16).attr("y",9)
            .attr("fill", hidden ? "#999" : "#333")
            .style("text-decoration", hidden ? "line-through" : "none")
            .text(g.name);
      }});
    }}

    function draw() {{
      const xKey = xSelect.property("value");
      const yKey = ySelect.property("value");
      const xMeta = metrics.find(m => m.key === xKey) || {{}};
      const yMeta = metrics.find(m => m.key === yKey) || {{}};

      document.getElementById("chart-title").textContent =
        `${{yMeta.label || yKey}} vs ${{xMeta.label || xKey}}`;

      const data = records.filter(r =>
        r.metrics[xKey] !== null && r.metrics[xKey] !== undefined &&
        r.metrics[yKey] !== null && r.metrics[yKey] !== undefined &&
        !hiddenGroups.has(groupNameOf(r)));

      const xExtent = d3.extent(data, d => d.metrics[xKey]);
      const yExtent = d3.extent(data, d => d.metrics[yKey]);
      const xPad = (xExtent[1] - xExtent[0]) * 0.1 || 1;
      const yPad = (yExtent[1] - yExtent[0]) * 0.1 || 1;
      x.domain([xExtent[0]-xPad, xExtent[1]+xPad]);
      y.domain([yExtent[0]-yPad, yExtent[1]+yPad]);

      gGrid.call(d3.axisLeft(y).tickSize(-W).tickFormat(""));
      gXAxis.call(d3.axisBottom(x).ticks(8).tickFormat(d3.format(",.2~f")));
      gYAxis.call(d3.axisLeft(y).ticks(8).tickFormat(d3.format(",.2~f")));
      xLabel.text(xMeta.label || xKey);
      yLabel.text(yMeta.label || yKey);

      dotsG.selectAll("circle").data(data, d => d.label)
        .join("circle")
        .attr("class", "dot")
        .attr("r", 6)
        .attr("fill", d => d.colors[colorMode])
        .attr("cx", d => x(d.metrics[xKey]))
        .attr("cy", d => y(d.metrics[yKey]))
        .on("mousemove", function(event, d) {{
          const xd = xMeta.decimals ?? 2, yd = yMeta.decimals ?? 2;
          tooltip.transition().duration(40).style("opacity",1);
          tooltip.html(
            `<strong>${{d.label}}</strong><br/>` +
            `Sample: ${{d.sample}}<br/>` +
            `Run: ${{d.run}}<br/>` +
            `Cell ID: ${{d.cellId}}<br/>` +
            `${{xMeta.label}}: <b>${{d3.format(",." + xd + "f")(d.metrics[xKey])}}</b><br/>` +
            `${{yMeta.label}}: <b>${{d3.format(",." + yd + "f")(d.metrics[yKey])}}</b>`
          ).style("left",(event.pageX+15)+"px")
           .style("top",Math.min(event.pageY-28, window.innerHeight-160)+"px");
        }})
        .on("mouseout", () => tooltip.transition().duration(200).style("opacity",0));
    }}

    xSelect.on("change", draw);
    ySelect.on("change", draw);
    colorSelect.on("change", function() {{
      colorMode = this.value;
      hiddenGroups = new Set();
      renderLegend();
      draw();
    }});

    if (metrics.length) {{
      xSelect.property("value", metrics.find(m => m.key === "loading_conc_pM") ? "loading_conc_pM" : metrics[0].key);
      ySelect.property("value", metrics.find(m => m.key === "hifi_yield_gb") ? "hifi_yield_gb" : (metrics[1] || metrics[0]).key);
      renderLegend();
      draw();
    }}
  </script>
</body>
</html>"""


def _html_escape_attr(s):
    return (s.replace("&", "&amp;")
             .replace('"', "&quot;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace("'", "&#x27;"))


# ---------------------------------------------------------------------------
# MAIN GENERATION FUNCTION
# ---------------------------------------------------------------------------

def generate_report(files, title, out_path):
    cells = load_all_cells(files)
    if not cells:
        sys.exit("ERROR: No cell rows found in the given file(s).")

    assign_colors(cells)

    n_runs = len({c["run_id"] or c["run_name"] for c in cells})
    print(f"Loaded {len(cells)} SMRT cell(s) across {n_runs} run(s) from {len(files)} file(s).")

    payload = build_chart_payload(cells)
    json_str = json.dumps(payload, separators=(",", ":"))

    bar_iframe = _html_escape_attr(BAR_CHART_IFRAME.format(plot_data_json=json_str))
    scatter_iframe = _html_escape_attr(SCATTER_IFRAME.format(plot_data_json=json_str))

    table_rows = "\n".join(render_table_row(rec) for rec in cells)
    table_footer = render_table_footer(cells)

    report = HTML_TEMPLATE.format(
        title=title,
        n_cells=len(cells),
        n_runs=n_runs,
        table_header=render_table_header(),
        table_rows=table_rows,
        table_footer=table_footer,
        filter_run_options=render_filter_options(cells, "run_name"),
        filter_cell_options=render_filter_options(cells, "cell_id"),
        filter_sample_options=render_filter_options(cells, "sample_name"),
        bar_iframe=bar_iframe,
        scatter_iframe=scatter_iframe,
        generated_on=datetime.now().strftime("%Y-%m-%d %H:%M"),
    )

    Path(out_path).write_text(report, encoding="utf-8")
    print(f"\nReport written to: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive PacBio HiFi run QC HTML report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="run-qc-export.csv file(s). Auto-discovered if omitted.")
    parser.add_argument("--title", default="HiFi Sequencing Report",
                        help="Report title (default: 'HiFi Sequencing Report')")
    parser.add_argument("--out", default="report.html",
                        help="Output HTML file (default: report.html)")
    args = parser.parse_args()

    files = sorted(args.files or glob.glob("*qc-export*.csv"))
    if not files:
        sys.exit("ERROR: No input files found. Pass file(s) explicitly or place "
                 "*qc-export*.csv files in the current directory.")

    generate_report(files=files, title=args.title, out_path=args.out)


if __name__ == "__main__":
    main()
