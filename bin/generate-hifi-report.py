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
  --out     PATH   Output HTML file (default: hifi-report.html).
  --demux   FILE.. SMRT Link "*.barcodes.summary.csv" demultiplexing file(s);
                   adds a Demultiplexing section. Such files can also be
                   passed as FILE arguments (they are recognised by their
                   columns).

If no FILE arguments are given, the script auto-discovers files matching
"*qc-export*.csv" and "*barcodes.summary*.csv" in the current directory.

Each barcodes summary is matched to one acquisition: its file name
("analysis-<sample>-<job>-m<instrument>_<YYMMDD>_<HHMMSS>_s<N>...") gives the
instrument, movie date and sample name, which must agree with the run-qc-export
row; when several rows still fit, the one whose HiFi yield is closest to the
barcodes' total yield wins.

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
import math
import re
import sys
from datetime import datetime
from pathlib import Path

from report_common import (D3_SOURCE_TAG, PAGE_TIP, TABLE_CSV_JS, fill, footer_html,
                           load_d3, render_chart_iframe, tip_attr, tip_html)

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
        f = float(v.replace(",", ""))
    except ValueError:
        return None
    # "nan"/"inf" parse as floats but mean the value is missing.
    return f if math.isfinite(f) else None


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
# by key, in display order). "read_lengths" and "loading" are combined mini
# charts (see READ_LENGTH_MARKS / LOADING_PARTS). All metrics remain
# available in the plot dropdowns regardless of this list.
TABLE_METRIC_KEYS = [
    "hifi_yield_gb", "q20_reads", "read_lengths", "read_quality_q", "q30_bq_pct", "loading",
]

# Shorter table headers for some metrics (CSV export keeps the full label).
SHORT_LABELS = {"read_quality_q": "Median Q", "q30_bq_pct": "Q30+ BQ"}

# Zero-based amounts that get a data bar in the table.
BAR_METRIC_KEYS = {"hifi_yield_gb", "q20_reads"}

# Read-length mini chart: (metric key, label, marker colour), on a kb axis
# shared by all rows. Insert size is drawn as a ring and mean length as a dot
# on top, since the two usually overlap.
READ_LENGTH_MARKS = [
    ("polymerase_rl_bp",  "Polymerase RL", "#64748b"),
    ("polymerase_n50_bp", "Polymerase N50", "#1d4ed8"),
    ("insert_size_bp",    "Insert size",   "#d97706"),
    ("mean_length_bp",    "Mean length",   "#059669"),
]
RING_MARKS = {"insert_size_bp"}
# Dashed guide lines through the read-length column, in kb.
READ_LENGTH_GUIDES_KB = [10, 100, 200]
# Read-length mini chart geometry (px); the guides are positioned from it.
RL_CHART_W, RL_CHART_PAD, RL_CELL_PAD = 124, 5, 10

# Loading mini chart: ZMW occupancy as a 100% stacked bar.
LOADING_PARTS = [
    ("p0_pct", "P0 (empty)",    "#cbd5e0"),
    ("p1_pct", "P1 (single)",   "#22c55e"),
    ("p2_pct", "P2 (multiple)", "#f59e0b"),
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

    # Labels ("<run date> <well>") name an acquisition in the charts and must be
    # unique: when two runs share a date, use the full run name instead.
    short = {}
    for rec in cells:
        short.setdefault(rec["label"], set()).add(rec["run_name"])
    for rec in cells:
        if len(short[rec["label"]]) > 1:
            rec["label"] = f"{rec['run_name']} {rec['well']}".strip()
    return cells


# ---------------------------------------------------------------------------
# DEMULTIPLEXING (SMRT Link "*.barcodes.summary.csv", one row per barcode)
# ---------------------------------------------------------------------------

# SMRT Link names these files "analysis-<sample>-<job id>-<movie>.barcodes.summary.csv",
# the movie being "m<instrument>_<YYMMDD>_<HHMMSS>_s<N>".
DEMUX_NAME_RE = re.compile(
    r"^(?:analysis-)?(?P<sample>.+?)-(?P<job>\d+)-"
    r"(?P<movie>m(?P<inst>\w+?)_(?P<date>\d{6})_(?P<time>\d{6})_s\d+)\.barcodes\.summary", re.I)
# Run IDs look like "r84180_20261021_131658".
RUN_ID_RE = re.compile(r"^r\w+?_(\d{8})_(\d{6})$")
# A movie starts after its run does; allow time-zone slack before and a run
# with several cells in sequence after.
MOVIE_WINDOW_DAYS = (-1, 3)
# Warn when the barcodes' summed yield is this far off the matched cell's.
YIELD_MISMATCH = 0.05
# Barcode rows that are not a sample (reads with no barcode found).
NOT_BARCODED = "Not Barcoded"


def is_demux_csv(path):
    """True for a barcodes summary CSV (as opposed to a run-qc-export)."""
    try:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            header = next(csv.reader(fh), [])
    except OSError:
        return False
    return "Barcode" in header and "HiFi Reads" in header


def _norm_name(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def load_demux(path):
    name = Path(path).name
    rows = []
    for raw in parse_qc_csv(path):
        barcode = (raw.get("Barcode") or "").strip()
        yield_bp = _float(raw.get("HiFi Yield (bp)"))
        rows.append({
            "bio_sample": (raw.get("Sample Name") or "").strip(),
            "barcode":    barcode,
            "barcoded":   barcode.lower() != NOT_BARCODED.lower(),
            "bc_quality": _float(raw.get("Barcode Quality")),
            "reads":      _int(raw.get("HiFi Reads")),
            "mean_len":   _int(raw.get("HiFi Read Length (mean, bp)")),
            "qv":         _q_value(raw.get("HiFi Read Quality (median, QV)")),
            "yield_gb":   round(yield_bp / 1e9, 3) if yield_bp is not None else None,
        })
    # Barcoded samples first (in file order), then the unbarcoded reads.
    rows.sort(key=lambda r: not r["barcoded"])
    total_reads = sum(r["reads"] or 0 for r in rows)
    for r in rows:
        r["reads_pct"] = round(100 * (r["reads"] or 0) / total_reads, 2) if total_reads else None

    m = DEMUX_NAME_RE.match(name)
    info = m.groupdict() if m else {}
    movie_dt = None
    if m:
        try:
            movie_dt = datetime.strptime(info["date"] + info["time"], "%y%m%d%H%M%S")
        except ValueError:
            pass
    return {"file": name, "rows": rows, "cell": None,
            "movie": info.get("movie", ""), "instrument": info.get("inst", ""),
            "sample_hint": info.get("sample", ""), "movie_dt": movie_dt,
            "yield_gb": sum(r["yield_gb"] or 0 for r in rows)}


def _run_dt(rec):
    """Run start from the Run ID (instrument clock, like the movie name)."""
    m = RUN_ID_RE.match(rec["run_id"])
    if m:
        try:
            return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")
        except ValueError:
            pass
    return rec["run_start_dt"]


def match_demux(demux, cells):
    """Pick the acquisition a barcodes summary belongs to.

    Instrument and movie date (from the file name) must agree with the run;
    among those, cells whose sample name is in the file name are preferred,
    and the cell with the closest HiFi yield to the barcodes' total wins.
    Without a parsable file name all cells are candidates.
    """
    cands = list(cells)
    if demux["instrument"]:
        cands = [c for c in cands if not c["instrument_sn"] or c["instrument_sn"] == demux["instrument"]]
    if demux["movie_dt"]:
        lo, hi = MOVIE_WINDOW_DAYS
        def in_window(c):
            dt = _run_dt(c)
            return dt is None or lo <= (demux["movie_dt"] - dt).total_seconds() / 86400 <= hi
        cands = [c for c in cands if in_window(c)]
    hint = _norm_name(demux["sample_hint"] or demux["file"])
    named = [c for c in cands if c["sample_name"] and _norm_name(c["sample_name"]) in hint]
    cands = named or cands
    if not cands:
        return None

    def yield_gap(c):
        y = c["metrics"].get("hifi_yield_gb")
        return abs(y - demux["yield_gb"]) / max(y, 1e-9) if y else float("inf")
    best = min(cands, key=yield_gap)
    # With nothing but yield to go on, a poor yield match is no match.
    if not named and not (demux["instrument"] or demux["movie_dt"]) and yield_gap(best) > YIELD_MISMATCH:
        return None
    if yield_gap(best) > YIELD_MISMATCH:
        print(f"WARNING: {demux['file']}: barcode yield {demux['yield_gb']:.2f} Gb differs from "
              f"matched cell yield {best['metrics'].get('hifi_yield_gb')} Gb.")
    return best


def load_all_demux(paths, cells):
    demuxes, seen = [], set()
    for path in paths:
        d = load_demux(path)
        key = d["movie"] or d["file"]
        if key in seen:
            continue
        seen.add(key)
        if not d["rows"]:
            print(f"WARNING: {d['file']}: no barcode rows, skipped.")
            continue
        d["cell"] = match_demux(d, cells)
        if d["cell"]:
            print(f"Demux {d['file']} -> {d['cell']['run_name']} {d['cell']['well']} "
                  f"({d['cell']['sample_name']}, SMRT Cell {d['cell']['cell_id']})")
        else:
            print(f"WARNING: {d['file']}: no matching acquisition found; shown as unmatched.")
        demuxes.append(d)
    taken = {}
    for d in demuxes:
        if d["cell"]:
            taken.setdefault(id(d["cell"]), []).append(d["file"])
    for files in taken.values():
        if len(files) > 1:
            print(f"WARNING: several barcode files matched the same acquisition: {', '.join(files)}")
    order = {id(c): i for i, c in enumerate(cells)}
    demuxes.sort(key=lambda d: order.get(id(d["cell"]), len(cells)))
    return demuxes


# Group name for cells without a sample name (never the run name, which
# would mix runs into the sample grouping).
NO_SAMPLE = "(no sample name)"

# Color-by grouping modes, in dropdown display order. Key is used internally
# / in JS; label is the human-readable dropdown option text.
COLOR_GROUP_MODES = [
    ("sample", "Sample"),
    ("run", "Run"),
    ("cell", "SMRT Cell"),
]


def _group_value(rec, mode):
    """Value used to group/color a cell for a given color-by mode."""
    if mode == "sample":
        return rec["sample_name"] or NO_SAMPLE
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
# dropdown), metric columns (numeric, pulled from rec["metrics"]), the derived
# "Acquisition" column (Use + Is Multi-Use Cell combined) and the two mini
# chart columns.
COLUMNS = [
    {"key": "run_name",    "label": "Run Name", "type": "text", "filter": True},
    {"key": "cell_id",     "label": "SMRT Cell", "type": "text", "filter": True, "csv_cols": "SMRT Cell ID"},
    {"key": "acquisition", "label": "Use", "type": "number", "csv_cols": "Acquisition",
     "dec": 0, "filter": False, "title": "Which use of the SMRT Cell this row is (1 = first use); "
                                          "each row is one acquisition, i.e. one use of a SMRT Cell in a well of a run."},
    {"key": "sample_name", "label": "Sample",    "type": "text", "filter": True},
    {"key": "well",        "label": "Well",      "type": "text", "filter": False},
]
for _mkey in TABLE_METRIC_KEYS:
    if _mkey == "read_lengths":
        COLUMNS.append({
            "key": _mkey, "label": "Read Lengths (kb)", "type": "number", "filter": False,
            "csv_cols": "|".join(METRICS_BY_KEY[k][0] for k, _l, _c in READ_LENGTH_MARKS),
            "title": "Markers on a kb axis shared by all rows: "
                     + ", ".join(label for _k, label, _c in READ_LENGTH_MARKS)
                     + ". Dashed lines at " + ", ".join(f"{kb}" for kb in READ_LENGTH_GUIDES_KB)
                     + " kb. Sorts by polymerase N50."})
    elif _mkey == "loading":
        COLUMNS.append({
            "key": _mkey, "label": "Loading", "type": "number", "filter": False,
            "csv_cols": "|".join(METRICS_BY_KEY[k][0] for k, _l, _c in LOADING_PARTS),
            "title": "Share of ZMWs that are empty (P0), hold one polymerase (P1) "
                     "or several (P2). Sorts by P1."})
    else:
        _label, _dec = METRICS_BY_KEY[_mkey]
        _col = {"key": _mkey, "label": SHORT_LABELS.get(_mkey, _label), "type": "number",
                "dec": _dec, "bar": _mkey in BAR_METRIC_KEYS, "filter": False}
        if _mkey in SHORT_LABELS:
            _col["csv_cols"] = _label
        COLUMNS.append(_col)

FILTER_COLUMNS = [c for c in COLUMNS if c["filter"]]
# Columns summed in the run subtotal and Total rows.
SUM_KEYS = ["hifi_yield_gb", "q20_reads"]


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
    if key == "read_lengths":
        return rec["metrics"].get("polymerase_n50_bp")
    if key == "loading":
        return rec["metrics"].get("p1_pct")
    return rec["metrics"].get(key)


def _col_display(rec, col):
    val = _col_value(rec, col)
    if col["type"] == "text":
        return val or ""
    if col["key"] == "acquisition":
        if val is None:
            return "-"
        return f"{val:,.0f} (multi-use)" if rec["multi_use"] else f"{val:,.0f}"
    text = _fmt_number(col["dec"], val)
    # Percent metrics whose short header dropped the "(%)" carry the sign instead.
    return text + "%" if col["key"] in SHORT_LABELS and col["key"].endswith("_pct") and val is not None else text


DOT = '<span style="display:inline-block;width:11px;height:11px;border-radius:50%;background:#9ca3af;margin:0 2px;vertical-align:middle;"></span>'

def _acq_dots(val):
    """Render acquisition number as 1–N big grey dots (raw HTML)."""
    if val is None:
        return "-"
    n = max(1, min(int(val), 8))  # cap at 8 for safety
    return DOT * n


def _rl_px(bp, max_bp):
    """x of a length (bp) in the read-length chart; axis shared by all rows."""
    return RL_CHART_PAD + (RL_CHART_W - 2 * RL_CHART_PAD) * bp / max(max_bp, 1)


def read_length_guide_css(max_bp):
    """CSS drawing the kb guides as dashed lines filling each cell row.

    They are a background of the column's cells in the sample rows, so each
    line spans the full row height. The chart is right-aligned in the cell,
    so each line is placed from the cell's right edge.
    """
    offsets = [RL_CELL_PAD + RL_CHART_W - _rl_px(kb * 1000, max_bp) - 0.5
               for kb in READ_LENGTH_GUIDES_KB if kb * 1000 <= max_bp]
    if not offsets:
        return ""
    dash = "linear-gradient(#a0aec0 50%, transparent 50%)"
    return ("#cell-table tr.cell-row td.rl-col{"
            f"background-image:{', '.join([dash] * len(offsets))};"
            f"background-size:1px 4px;background-repeat:repeat-y;"
            f"background-position:{', '.join(f'right {o:.1f}px top 0' for o in offsets)};}}")


def _max_read_bp(cells):
    return max((rec["metrics"].get(k) or 0 for rec in cells
                for k, _l, _c in READ_LENGTH_MARKS), default=1)


def _read_length_chart(rec, max_bp, h=22):
    """Read-length markers on a kb axis shared by all rows (inline SVG)."""
    w, pad = RL_CHART_W, RL_CHART_PAD
    def px(bp): return _rl_px(bp, max_bp)
    marks = []
    for key, label, color in READ_LENGTH_MARKS:
        v = rec["metrics"].get(key)
        if v is None:
            continue
        if key in RING_MARKS:
            marks.append(f'<circle cx="{px(v):.1f}" cy="{h / 2}" r="5.5" fill="#fff" '
                         f'stroke="{color}" stroke-width="2"/>')
        else:
            marks.append(f'<circle cx="{px(v):.1f}" cy="{h / 2}" r="3.5" fill="{color}" '
                         f'stroke="#fff" stroke-width="1"/>')
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" '
            f'aria-label="Read lengths">'
            f'<line x1="{pad}" x2="{w - pad}" y1="{h / 2}" y2="{h / 2}" stroke="#e2e8f0" '
            f'stroke-width="2"/>{"".join(marks)}</svg>')


def _loading_bar(rec, w=64, h=12):
    """P0/P1/P2 as a 100% stacked bar, with P1 printed next to it."""
    parts = [(label, rec["metrics"].get(key), color) for key, label, color in LOADING_PARTS]
    if all(v is None for _l, v, _c in parts):
        return "-"
    total = sum(v or 0 for _l, v, _c in parts) or 1
    x, rects = 0.0, []
    for _label, v, color in parts:
        width = w * (v or 0) / total
        rects.append(f'<rect x="{x:.2f}" y="0" width="{width:.2f}" height="{h}" fill="{color}"/>')
        x += width
    label_text = ", ".join(f"{label.split()[0]} {v:.1f}%" for label, v, _c in parts if v is not None)
    p1 = rec["metrics"].get("p1_pct")
    p1_text = "" if p1 is None else f'<span class="load-p1">{p1:.1f}%</span>'
    return (f'<span class="load"><svg width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img" aria-label="Loading {label_text}">'
            f'{"".join(rects)}</svg>{p1_text}</span>')


def _read_length_tip(rec):
    # Ring markers get an outlined swatch, like in the chart.
    rows = [(f"ring:{c}" if k in RING_MARKS else c, label, f"{rec['metrics'][k] / 1000:,.1f} kb")
            for k, label, c in READ_LENGTH_MARKS if rec["metrics"].get(k) is not None]
    return tip_html(f"{_html_escape(rec['label'])}: read lengths", rows)


def _loading_tip(rec):
    rows = [(c, label, f"{rec['metrics'][k]:.1f}%")
            for k, label, c in LOADING_PARTS if rec["metrics"].get(k) is not None]
    return tip_html(f"{_html_escape(rec['label'])}: loading", rows,
                    "Share of ZMWs that are empty, hold one polymerase, or several")


def _csv_values(rec, keys):
    vals = [rec["metrics"].get(k) for k in keys]
    return "|".join("" if v is None else repr(v) for v in vals)


def render_table_header():
    ths = []
    for col in COLUMNS:
        attrs = [f'data-type="{col["type"]}"', f'data-key="{col["key"]}"', 'onclick="sortTable(this)"']
        if col["type"] == "number":
            attrs.append('style="text-align:right"')
        if col.get("bar"):
            attrs.append('data-bar="true"')
        if col.get("csv_cols"):
            attrs.append(f'data-csv-cols="{_html_escape(col["csv_cols"])}"')
        if col.get("title"):
            attrs.append(tip_attr(tip_html(_html_escape(col["label"]), text=_html_escape(col["title"]))))
        ths.append(f'<th {" ".join(attrs)}>{col["label"]}<span class="sort-arrow"></span></th>')
    return "".join(ths)


# Text columns holding a "Color by" group, by mode.
TAG_KEYS = {"run_name": "run", "cell_id": "cell", "sample_name": "sample"}


def render_table_row(rec, max_read_bp):
    tds = []
    for i, col in enumerate(COLUMNS):
        val = _col_value(rec, col)
        if col["type"] == "number":
            data_value = "" if val is None else repr(val)
        else:
            data_value = (val or "").lower()

        style_parts = []
        class_parts = []
        extra = ""
        if col["type"] == "number":
            style_parts.append("text-align:right")
        if i == 0:
            class_parts.append("color-stripe")
            # An inset shadow, not a border, so the stripe moves with the sticky cell.
            style_parts.append(f"box-shadow:inset 4px 0 0 {rec['colors']['sample']},inset -1px 0 0 #e2e8f0")

        if col["key"] == "acquisition":
            cell_content = _acq_dots(val)
            if val is not None:
                extra = " " + tip_attr(tip_html(f"Use {val}", text=(
                    f"Use {val} of SMRT Cell {_html_escape(rec['cell_id'])}"
                    + (" (multi-use SMRT Cell)" if rec["multi_use"] else ""))))
        elif col["key"] == "read_lengths":
            cell_content = _read_length_chart(rec, max_read_bp)
            class_parts += ["mini", "rl-col"]
            extra = (f' data-csv="{_csv_values(rec, [k for k, _l, _c in READ_LENGTH_MARKS])}" '
                     + tip_attr(_read_length_tip(rec)))
        elif col["key"] == "loading":
            cell_content = _loading_bar(rec)
            class_parts.append("mini")
            extra = (f' data-csv="{_csv_values(rec, [k for k, _l, _c in LOADING_PARTS])}" '
                     + tip_attr(_loading_tip(rec)))
        elif col["key"] in TAG_KEYS:
            # Tagged so the page can colour the column chosen in "Color by".
            cell_content = f'<span class="tag">{_html_escape(_col_display(rec, col))}</span>'
        else:
            cell_content = _html_escape(_col_display(rec, col))
        style_attr = f' style="{";".join(style_parts)}"' if style_parts else ""
        class_attr = f' class="{" ".join(class_parts)}"' if class_parts else ""
        tds.append(f'<td{class_attr}{style_attr} data-value="{_html_escape(data_value)}"{extra}>'
                   f'{cell_content}</td>')

    row_attrs = (f'class="cell-row" '
                f'data-run="{_html_escape(rec["run_name"])}" '
                f'data-cell="{_html_escape(rec["cell_id"])}" '
                f'data-sample="{_html_escape(rec["sample_name"])}" '
                f'data-color-sample="{rec["colors"]["sample"]}" '
                f'data-color-run="{rec["colors"]["run"]}" '
                f'data-color-cell="{rec["colors"]["cell"]}"')
    return f"    <tr {row_attrs}>" + "".join(tds) + "</tr>"


def _sum_row(cls, label, csv_label):
    """Row with label + summed columns; values are filled in by the page JS."""
    tds = []
    for i, col in enumerate(COLUMNS):
        if i == 0:
            tds.append(f'<td class="sum-label" data-csv="{csv_label}">{label}</td>')
        elif col["key"] in SUM_KEYS:
            tds.append(f'<td style="text-align:right" data-sum="{col["key"]}" data-value=""></td>')
        else:
            tds.append("<td></td>")
    return f'    <tr class="{cls}">' + "".join(tds) + "</tr>"


def render_table_body(cells):
    """One <tbody> per run (cells in run order) ending in a subtotal row."""
    max_read_bp = _max_read_bp(cells)
    runs = {}
    for rec in cells:
        runs.setdefault(rec["run_name"] or "Unknown", []).append(rec)
    bodies = []
    for run, recs in runs.items():
        rows = "\n".join(render_table_row(rec, max_read_bp) for rec in recs)
        bodies.append(f'  <tbody class="run-group" data-run="{_html_escape(run)}">\n{rows}\n'
                      + _sum_row("subtotal", "Subtotal", "Subtotal") + "\n  </tbody>")
    return "\n".join(bodies)


def render_filter_options(cells, key):
    """Checkbox list for one multi-select filter dropdown."""
    values = sorted({rec[key] for rec in cells if rec[key]})
    return "".join(f'<label><input type="checkbox" value="{_html_escape(v)}">{_html_escape(v)}</label>'
                   for v in values)


def render_table_footer():
    return _sum_row("total-row", "<strong>Total</strong>", "Total")


def render_table_legend():
    """Key for the two mini chart columns, shown under the table."""
    def swatch(color, shape):
        radius = "2px" if shape == "box" else "50%"
        fill = (f"background:#fff;border:2px solid {color};width:7px;height:7px;"
                if shape == "ring" else f"background:{color};width:9px;height:9px;")
        return (f'<span style="display:inline-block;{fill}border-radius:{radius};'
                f'margin:0 4px 0 10px;vertical-align:middle;"></span>')
    lengths = "".join(swatch(c, "ring" if k in RING_MARKS else "dot") + label
                      for k, label, c in READ_LENGTH_MARKS)
    loading = "".join(swatch(c, "box") + label for _k, label, c in LOADING_PARTS)
    guides = ", ".join(str(kb) for kb in READ_LENGTH_GUIDES_KB)
    return (f"<b>Read lengths:</b>{lengths} &nbsp;(dashed lines: {guides} kb) "
            f"&nbsp;&nbsp; <b>Loading:</b>{loading}")


# ---------------------------------------------------------------------------
# DEMUX TABLE: one <tbody> per acquisition, one row per barcode
# ---------------------------------------------------------------------------

# (key, label, decimals, data bar, CSV column(s)); text columns have decimals
# None. Short labels keep the table narrow; the CSV gets the full names.
DEMUX_COLUMNS = [
    ("acq",        "Acquisition",  None, False, "Run Name|Well|SMRT Cell"),
    ("bio_sample", "Bio Sample",   None, False, None),
    ("barcode",    "Barcode",      None, False, None),
    ("bc_quality", "BC Quality",   1,    False, "Barcode Quality"),
    ("reads",      "HiFi Reads",   0,    True,  None),
    ("reads_pct",  "Reads %",      1,    False, "Reads (%)"),
    ("mean_len",   "Mean Length",  0,    False, "Mean Length (bp)"),
    ("qv",         "Median QV",    0,    False, None),
    ("yield_gb",   "Yield (Gb)",   2,    True,  "HiFi Yield (Gb)"),
]
DEMUX_SUM_KEYS = {"reads", "reads_pct", "yield_gb"}


def _demux_ident(d):
    c = d["cell"]
    if c is None:
        return {"run_name": "(unmatched)", "well": "", "cell_id": "", "sample": ""}
    return {"run_name": c["run_name"], "well": c["well"], "cell_id": c["cell_id"],
            "sample": c["sample_name"]}


def _demux_label(d):
    return d["cell"]["label"] if d["cell"] else d["file"]


def barcode_cv(rows):
    """Coefficient of variation (%) of reads across barcoded samples."""
    vals = [r["reads"] for r in rows if r["barcoded"] and r["reads"] is not None]
    if len(vals) < 2:
        return None
    mean = sum(vals) / len(vals)
    sd = (sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
    return 100 * sd / mean if mean else None


def render_demux_header():
    ths = []
    for key, label, dec, _bar, csv_cols in DEMUX_COLUMNS:
        num = dec is not None
        attrs = [f'data-type="{"number" if num else "text"}"']
        if num:
            attrs.append('style="text-align:right"')
        if csv_cols:
            attrs.append(f'data-csv-cols="{_html_escape(csv_cols)}"')
        if key == "reads_pct":
            attrs.append(tip_attr(tip_html("Reads %", text="Share of the acquisition's HiFi reads.")))
        elif key == "bc_quality":
            attrs.append(tip_attr(tip_html("Barcode Quality", text="Mean barcode score (0-100).")))
        ths.append(f'<th {" ".join(attrs)}>{label}</th>')
    return "".join(ths)


# "Color by" in the demux table and chart: the acquisition's groups, as in
# the Acquisition Summary.
DEMUX_COLOR_MODES = COLOR_GROUP_MODES
UNMATCHED_COLOR = "#cbd5e0"


def demux_colors(demuxes):
    """Set r["colors"] / r["names"] {mode: ...} on every barcode row."""
    for d in demuxes:
        for r in d["rows"]:
            r["colors"] = {mode: d["cell"]["colors"][mode] if d["cell"] else UNMATCHED_COLOR
                           for mode, _label in DEMUX_COLOR_MODES}
            r["names"] = {mode: _group_value(d["cell"], mode) if d["cell"] else "(unmatched)"
                          for mode, _label in DEMUX_COLOR_MODES}


def render_demux_bodies(demuxes):
    bodies = []
    for d in demuxes:
        # Data bars are scaled per acquisition, to show how evenly it split.
        maxes = {key: max((r[key] or 0 for r in d["rows"]), default=0)
                 for key, _l, _d, bar, _c in DEMUX_COLUMNS if bar}
        ident = _demux_ident(d)
        rows = []
        acq_csv = _html_escape("|".join(ident[k] for k in ("run_name", "well", "cell_id")))
        acq_tip = tip_attr(tip_html(_html_escape(_demux_label(d)), [
            (None, "Run", _html_escape(ident["run_name"])),
            (None, "Well", _html_escape(ident["well"] or "-")),
            (None, "SMRT Cell", _html_escape(ident["cell_id"] or "-")),
            (None, "Sample", _html_escape(ident["sample"] or "-")),
        ], _html_escape(d["file"])))
        for r in d["rows"]:
            tds = []
            for i, (key, _label, dec, bar, _csv) in enumerate(DEMUX_COLUMNS):
                val = r.get(key)
                styles = []
                cls = extra = ""
                if key == "acq":
                    cls = ' class="color-stripe"'
                    styles.append(f"box-shadow:inset 4px 0 0 {r['colors']['sample']},inset -1px 0 0 #e2e8f0")
                    extra = f' data-csv="{acq_csv}" {acq_tip}'
                    # Tagged parts are coloured by the matching "Color by" mode.
                    sub = " &middot; ".join(f'<span class="tag" data-tag="{mode}">{_html_escape(v)}</span>'
                                            for mode, v in (("sample", ident["sample"]), ("cell", ident["cell_id"]))
                                            if v)
                    content = (f'<span class="tag" data-tag="run">{_html_escape(_demux_label(d))}</span>'
                               f'<div class="acq-sub">{sub}</div>')
                    style = f' style="{";".join(styles)}"'
                    tds.append(f'<td{cls}{style}{extra}>{content}</td>')
                    continue
                if dec is None:
                    content, data_value = _html_escape(val or ""), _html_escape(val or "")
                else:
                    styles.append("text-align:right")
                    content = _fmt_number(dec, val)
                    if key == "qv" and val is not None:
                        content = f"Q{content}"
                    data_value = "" if val is None else repr(val)
                    if bar and val is not None and maxes[key]:
                        content += (f'<div class="bar-track"><div class="bar-fill" '
                                    f'style="width:{round(100 * val / maxes[key])}%"></div></div>')
                style = f' style="{";".join(styles)}"' if styles else ""
                tds.append(f'<td{cls}{style} data-value="{data_value}"{extra}>{content}</td>')
            row_cls = "demux-row" + ("" if r["barcoded"] else " unbarcoded")
            mode_attrs = "".join(f' data-color-{m}="{r["colors"][m]}" data-name-{m}="{_html_escape(r["names"][m])}"'
                                 for m, _l in DEMUX_COLOR_MODES)
            rows.append(f'    <tr class="{row_cls}"{mode_attrs}>' + "".join(tds) + "</tr>")

        n_bc = sum(r["barcoded"] for r in d["rows"])
        cv = barcode_cv(d["rows"])
        cv_text = "" if cv is None else f" &middot; CV {cv:.0f}%"
        tds = []
        for i, (key, _label, dec, _bar, _csv) in enumerate(DEMUX_COLUMNS):
            if i == 0:
                tds.append(f'<td class="sum-label" ' + tip_attr(tip_html(
                    "Subtotal", text="CV: coefficient of variation of HiFi reads across the "
                                     "barcoded samples (lower is more even).")) +
                           f'>Subtotal <span class="sub">{n_bc} barcode{"" if n_bc == 1 else "s"}{cv_text}</span></td>')
            elif key in DEMUX_SUM_KEYS:
                total = sum(r[key] or 0 for r in d["rows"])
                tds.append(f'<td style="text-align:right">{_fmt_number(dec, total)}</td>')
            else:
                tds.append("<td></td>")
        rows.append('    <tr class="subtotal">' + "".join(tds) + "</tr>")
        # Identity of the acquisition, for the Acquisition Summary filters.
        ident_attrs = (f'data-run="{_html_escape(ident["run_name"])}" '
                       f'data-cell="{_html_escape(ident["cell_id"])}" '
                       f'data-sample="{_html_escape(ident["sample"])}"')
        bodies.append(f'  <tbody class="run-group" {ident_attrs}>\n' + "\n".join(rows) + "\n  </tbody>")
    return "\n".join(bodies)


def build_demux_payload(demuxes):
    groups = []
    for d in demuxes:
        ident = _demux_ident(d)
        groups.append({
            "label": _demux_label(d), "run": ident["run_name"], "well": ident["well"],
            "cellId": ident["cell_id"], "sample": ident["sample"],
            "rows": [{"name": r["bio_sample"] or r["barcode"], "barcode": r["barcode"],
                      "barcoded": r["barcoded"], "colors": r["colors"],
                      "metrics": {k: r[k] for k in ("reads", "reads_pct", "yield_gb",
                                                    "mean_len", "qv", "bc_quality")}}
                     for r in d["rows"]],
        })
    metrics = [
        {"key": "reads",      "label": "HiFi Reads",                  "decimals": 0, "even": True},
        {"key": "yield_gb",   "label": "HiFi Yield (Gb)",             "decimals": 2, "even": True},
        {"key": "reads_pct",  "label": "Reads (% of acquisition)",    "decimals": 1, "even": True},
        {"key": "mean_len",   "label": "HiFi Read Length (mean, bp)", "decimals": 0},
        {"key": "qv",         "label": "HiFi Read Quality (median, QV)", "decimals": 0},
        {"key": "bc_quality", "label": "Barcode Quality",             "decimals": 1},
    ]
    return {"groups": groups, "metrics": metrics}


def render_demux_section(demuxes):
    if not demuxes:
        return ""
    demux_colors(demuxes)
    n_rows = sum(len(d["rows"]) for d in demuxes)
    options = "".join(f'<option value="{m}">{label}</option>' for m, label in DEMUX_COLOR_MODES)
    iframe = _render_iframe(DEMUX_IFRAME, build_demux_payload(demuxes), "Demultiplexing")
    n_unmatched = sum(d["cell"] is None for d in demuxes)
    unmatched = (f" {n_unmatched} file{'' if n_unmatched == 1 else 's'} could not be matched to an "
                 "acquisition of the run-qc-export file(s) and "
                 f"{'is' if n_unmatched == 1 else 'are'} shown as \"(unmatched)\"." if n_unmatched else "")
    return f"""\
<div class="section-head">
  <h3>Demultiplexing</h3>
  <button class="btn" onclick="downloadTableCsv('demux-table', 'demux')">Download CSV</button>
</div>
<div class="filter-bar">
  <label for="demux-color-by">Color by</label>
  <select id="demux-color-by">{options}</select>
</div>
<div id="demux-legend" class="color-legend"></div>
<div class="table-wrap">
<table id="demux-table">
  <thead><tr>{render_demux_header()}</tr></thead>
{render_demux_bodies(demuxes)}
</table>
</div>
<p class="table-note">One row per barcode, grouped by acquisition; barcode files are matched to
acquisitions by instrument, movie date and sample name in the file name, then by HiFi yield.{unmatched}</p>
<div class="plot-container">
  <iframe id="demux-chart" srcdoc="{iframe}" height="{150 + 24 * n_rows + 60 * len(demuxes)}" frameborder="0" scrolling="no"></iframe>
</div>
<script>
  // Stripe, tinted tag and legend for the "Color by" mode; the chart follows.
  function applyDemuxColors(mode) {{
    const attr = 'color' + mode[0].toUpperCase() + mode.slice(1);
    const nameAttr = 'name' + mode[0].toUpperCase() + mode.slice(1);
    const seen = new Map();
    document.querySelectorAll('#demux-table tr.demux-row').forEach(r => {{
      const color = r.dataset[attr];
      r.querySelector('td.color-stripe').style.boxShadow = `inset 4px 0 0 ${{color}},inset -1px 0 0 #e2e8f0`;
      r.querySelectorAll('.tag').forEach(t => {{
        const on = t.dataset.tag === mode;
        t.classList.toggle('on', on);
        t.style.background = on ? `color-mix(in srgb, ${{color}} 16%, transparent)` : '';
        t.style.color = on ? `color-mix(in srgb, ${{color}} 75%, #000)` : '';
      }});
      if (!seen.has(r.dataset[nameAttr])) seen.set(r.dataset[nameAttr], color);
    }});
    const legend = document.getElementById('demux-legend');
    legend.innerHTML = '';
    seen.forEach((color, name) => {{
      const item = document.createElement('div');
      item.innerHTML = `<span style="background:${{color}}"></span>`;
      item.append(name);
      legend.appendChild(item);
    }});
    const chart = document.getElementById('demux-chart');
    try {{ chart.contentWindow.setColorMode && chart.contentWindow.setColorMode(mode); }} catch (e) {{}}
  }}
  const demuxColorBy = document.getElementById('demux-color-by');
  demuxColorBy.addEventListener('change', () => applyDemuxColors(demuxColorBy.value));
  applyDemuxColors(demuxColorBy.value);
  // The chart may load after this runs: give it the current mode once it is up.
  document.getElementById('demux-chart').addEventListener('load', () => applyDemuxColors(demuxColorBy.value));
</script>
"""


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
            "sample": rec["sample_name"] or NO_SAMPLE,
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
<title>%%TITLE%%</title>
<meta name="description" content="Interactive PacBio HiFi run QC report.">
%%D3_SOURCE%%
<style>
  *,*::before,*::after{box-sizing:border-box;}
  body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
        max-width:1160px;margin:0 auto;padding:24px 28px;color:#2d3748;
        line-height:1.55;background:#f7f8fa;}
  h2{color:#1a202c;border-bottom:2px solid #e2e8f0;padding-bottom:8px;margin-top:40px;}
  h3{color:#2d3748;margin-top:28px;font-size:1.05rem;}
  .subtitle{color:#718096;font-size:13.5px;margin-top:-8px;}
  .table-wrap{overflow-x:auto;border-radius:6px;box-shadow:0 1px 4px rgba(0,0,0,.07);}
  table{border-collapse:collapse;width:100%;margin:14px 0;font-size:13px;
         background:#fff;white-space:nowrap;}
  th{background:#edf2f7;color:#2d3748;font-weight:600;text-align:left;
      padding:9px 10px;border-bottom:2px solid #cbd5e0;position:sticky;top:0;
      cursor:pointer;user-select:none;}
  th:hover{background:#e2e8f0;}
  .sort-arrow{font-size:10px;color:#4a5568;margin-left:3px;}
  td{padding:8px 10px 14px;border-bottom:1px solid #e9ecef;position:relative;vertical-align:middle;}
  td.color-stripe{padding-left:14px;}
  tr:last-child td{border-bottom:none;}
  td.mini{padding:6px 10px;line-height:0;}
%%READ_LENGTH_GUIDE_CSS%%
  td.mini svg{vertical-align:middle;}
  .load{display:inline-flex;align-items:center;gap:6px;line-height:1;}
  .load svg{border-radius:3px;}
  .load-p1{font-size:12px;color:#4a5568;min-width:42px;text-align:right;}
  /* Run groups: subtotal row closes each run's <tbody> */
  tbody.run-group + tbody.run-group tr:first-child td{border-top:2px solid #e2e8f0;}
  tr.subtotal td{background-color:#f8fafc;font-weight:600;color:#4a5568;font-size:12.5px;
                 padding:6px 13px;border-bottom:none;}
  .sum-label .sub{font-weight:400;color:#718096;margin-left:4px;}
  /* Run Name column stays visible while the table scrolls sideways */
  #cell-table th:first-child,#cell-table td:first-child{position:sticky;left:0;z-index:2;
     box-shadow:inset -1px 0 0 #e2e8f0;}
  #cell-table tbody td:first-child{background:#fff;}
  #cell-table tr.subtotal td:first-child{background:#f8fafc;}
  .tag{display:inline-block;padding:1px 7px;margin:-1px -7px;border-radius:4px;}
  .tag.on{font-weight:600;}
  #demux-table th{cursor:default;}
  #demux-table th:hover{background:#edf2f7;}
  #demux-table th:first-child,#demux-table td:first-child{position:sticky;left:0;z-index:2;}
  #demux-table tbody td:first-child{background:#fff;}
  #demux-table tr.subtotal td:first-child{background:#f8fafc;}
  #demux-table tr.unbarcoded td{color:#94a3b8;}
  #demux-table .acq-sub .tag{padding:0 5px;margin:0 -5px;}
  .color-legend{display:flex;flex-wrap:wrap;gap:6px 12px;margin:2px 0 12px;font-size:12px;color:#555;}
  .color-legend div{display:flex;align-items:center;gap:5px;}
  .color-legend span{width:10px;height:10px;border-radius:3px;display:inline-block;}
  #demux-table .acq-sub{font-size:11.5px;color:#718096;line-height:1.3;}
  #demux-table tr.unbarcoded .bar-fill{background:#cbd5e0;}
  .table-note{font-size:12px;color:#718096;margin:6px 0 8px;}
  tr:hover td{background-color:rgba(59,130,246,.04);}
  .bar-track{position:absolute;bottom:4px;left:10px;right:10px;height:3px;
              background:rgba(59,130,246,.12);border-radius:99px;overflow:hidden;}
  .bar-fill{position:absolute;inset:0;height:100%;border-radius:99px;
             background:#3b82f6;
             transition:width .35s cubic-bezier(.4,0,.2,1);}
  tfoot tr td{background:#edf2f7;font-weight:600;border-top:2px solid #cbd5e0;
               color:#2d3748;padding:8px 11px;white-space:nowrap;}
  .filter-bar{display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin:14px 0 4px;}
  .filter-bar label{font-size:12.5px;font-weight:600;color:#555;}
  .filter-bar select{padding:5px 10px;font-size:12.5px;border:1px solid #cbd5e0;
                      border-radius:6px;background:#f8fafc;color:#2d3748;cursor:pointer;}
  /* Multi-select filter: a dropdown of checkboxes; none checked = all */
  .ms{position:relative;font-size:12.5px;}
  .ms summary{list-style:none;padding:5px 26px 5px 10px;border:1px solid #cbd5e0;border-radius:6px;
              background:#f8fafc;color:#2d3748;cursor:pointer;min-width:90px;max-width:220px;
              overflow:hidden;text-overflow:ellipsis;white-space:nowrap;position:relative;}
  .ms summary::-webkit-details-marker{display:none;}
  .ms summary::after{content:"\\25BE";position:absolute;right:9px;color:#718096;}
  .ms.active summary{border-color:#3b82f6;background:#eff6ff;}
  .ms-list{position:absolute;z-index:20;top:calc(100% + 4px);left:0;min-width:100%;max-height:300px;
           overflow-y:auto;background:#fff;border:1px solid #cbd5e0;border-radius:6px;
           box-shadow:0 4px 14px rgba(0,0,0,.12);padding:4px 0;}
  .ms-list label{display:flex;align-items:center;gap:7px;padding:4px 12px;white-space:nowrap;
                 font-weight:400;color:#2d3748;cursor:pointer;}
  .ms-list label:hover{background:#edf2f7;}
  .row-count{font-size:12.5px;color:#718096;margin-left:auto;}
  .plot-container{margin:20px 0;background:#fff;border-radius:8px;
                   box-shadow:0 1px 4px rgba(0,0,0,.08);border:1px solid #e2e8f0;overflow:hidden;}
  .plot-container iframe{border:none;display:block;width:100%;}
  .section-head{display:flex;align-items:center;justify-content:space-between;gap:12px;
                flex-wrap:wrap;margin-top:28px;}
  .section-head h3{margin:0;}
  .btn{font:inherit;font-size:12px;padding:4px 12px;border:1px solid #cbd5e0;border-radius:6px;
       background:#fff;color:#2d3748;cursor:pointer;}
  .btn:hover{background:#edf2f7;}
  .footer{margin-top:48px;border-top:1px solid #e2e8f0;padding-top:14px;
           color:#a0aec0;font-size:12px;text-align:center;}
  @media (max-width:640px){
    body{padding:16px;}
    h2{margin-top:16px;font-size:1.3rem;}
    .filter-bar{gap:8px 12px;}
    .row-count{margin-left:0;width:100%;}
  }
</style>
</head>
<body>

<h2>%%TITLE%%</h2>
<p class="subtitle">%%COUNTS%% &mdash; generated on %%GENERATED_ON%%</p>

<h3>HiFi Yield and Reads per Sample</h3>
<div class="plot-container">
  <iframe srcdoc="%%SAMPLE_IFRAME%%" height="%%SAMPLE_IFRAME_H%%" frameborder="0" scrolling="no"></iframe>
</div>

<div class="section-head">
  <h3>Acquisition Summary</h3>
  <button class="btn" onclick="downloadTableCsv('cell-table', 'cells')" data-tip="<div class='tt-text'>Rows hidden by the filters are left out.</div>">Download CSV</button>
</div>
<div class="filter-bar">
  <label>Run Name</label>
  <details class="ms" id="filter-run" data-key="run"><summary>All</summary><div class="ms-list">%%FILTER_RUN_OPTIONS%%</div></details>
  <label>SMRT Cell</label>
  <details class="ms" id="filter-cell" data-key="cell"><summary>All</summary><div class="ms-list">%%FILTER_CELL_OPTIONS%%</div></details>
  <label>Sample</label>
  <details class="ms" id="filter-sample" data-key="sample"><summary>All</summary><div class="ms-list">%%FILTER_SAMPLE_OPTIONS%%</div></details>
  <button class="btn" id="filter-clear" type="button">Clear filters</button>
  <label for="color-by-select">Color by</label>
  <select id="color-by-select">
    <option value="sample">Sample</option>
    <option value="run">Run</option>
    <option value="cell">SMRT Cell</option>
  </select>
  <span class="row-count" id="row-count"></span>
</div>
<div id="table-legend" style="display:flex;flex-wrap:wrap;gap:12px;margin:2px 0 12px;font-size:12px;color:#555;"></div>
<div class="table-wrap">
<table id="cell-table">
  <thead><tr>%%TABLE_HEADER%%</tr></thead>
%%TABLE_BODIES%%
  <tfoot>
%%TABLE_FOOTER%%
  </tfoot>
</table>
</div>
<p class="table-note">%%TABLE_LEGEND%%<br>Each row is one acquisition: one use of a SMRT Cell, in a well of a run. Subtotals and the Total row count only the acquisitions shown by the filters. The filters also apply to every chart and to the Demultiplexing section.</p>
<script>
%%TABLE_CSV_JS%%

  const CELL_ROWS = 'tr.cell-row';
  const visibleCellRows = scope => Array.from(scope.querySelectorAll(CELL_ROWS))
    .filter(r => r.style.display !== 'none');
  const HEADER_THS = Array.from(document.querySelectorAll('#cell-table thead th'));

  // Sorting keeps the run groups: rows are sorted within each run's <tbody>,
  // and its subtotal row stays last.
  function sortTable(th) {
    const colIndex = HEADER_THS.indexOf(th);
    const type = th.dataset.type;
    const dir = th.dataset.dir === 'asc' ? 'desc' : 'asc';
    HEADER_THS.forEach(h => { h.dataset.dir = ''; h.querySelector('.sort-arrow').textContent = ''; });
    th.dataset.dir = dir;
    th.querySelector('.sort-arrow').textContent = dir === 'asc' ? '\u25B2' : '\u25BC';

    document.querySelectorAll('#cell-table tbody.run-group').forEach(tbody => {
      const rows = Array.from(tbody.querySelectorAll(CELL_ROWS));
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
      const subtotal = tbody.querySelector('tr.subtotal');
      rows.forEach(r => tbody.insertBefore(r, subtotal));
    });
    applyDataBars();
  }

  // Fill a subtotal/total row from the visible cell rows in `scope`.
  function fillSums(row, rows) {
    row.querySelectorAll('td[data-sum]').forEach(td => {
      const ci = Array.from(row.children).indexOf(td);
      const sum = rows.reduce((acc, r) => acc + (parseFloat(r.children[ci].dataset.value) || 0), 0);
      const isYield = td.dataset.sum === 'hifi_yield_gb';
      td.dataset.value = isYield ? sum.toFixed(2) : String(Math.round(sum));
      td.textContent = sum.toLocaleString(undefined, isYield
        ? {minimumFractionDigits: 2, maximumFractionDigits: 2} : {maximumFractionDigits: 0});
    });
  }

  function updateSums() {
    document.querySelectorAll('#cell-table tbody.run-group').forEach(tbody => {
      const rows = visibleCellRows(tbody);
      const subtotal = tbody.querySelector('tr.subtotal');
      // A subtotal of one acquisition repeats its row, so only show it for two or more.
      subtotal.style.display = rows.length > 1 ? '' : 'none';
      subtotal.querySelector('.sum-label').innerHTML =
        `Subtotal <span class="sub">${rows.length} acquisitions</span>`;
      fillSums(subtotal, rows);
    });
    const all = visibleCellRows(document.getElementById('cell-table'));
    const total = document.querySelector('#cell-table tr.total-row');
    fillSums(total, all);
    total.querySelector('.sum-label').innerHTML =
      `<strong>Total</strong> <span class="sub">${all.length} acquisition${all.length === 1 ? '' : 's'}</span>`;
  }

  // Filters: {run: [...], cell: [...], sample: [...]}; an empty list means all.
  // The chart iframes read it with parent.currentTableFilter() when they start
  // and are sent it with setFilter() whenever it changes.
  const FILTER_IDS = ['filter-run', 'filter-cell', 'filter-sample'];

  function currentTableFilter() {
    const f = {};
    FILTER_IDS.forEach(id => {
      const ms = document.getElementById(id);
      f[ms.dataset.key] = Array.from(ms.querySelectorAll('input:checked')).map(i => i.value);
    });
    return f;
  }
  const matchesFilter = (f, el) =>
    ['run', 'cell', 'sample'].every(k => !f[k].length || f[k].includes(el.dataset[k]));

  function updateFilterLabels() {
    FILTER_IDS.forEach(id => {
      const ms = document.getElementById(id);
      const checked = Array.from(ms.querySelectorAll('input:checked')).map(i => i.value);
      ms.querySelector('summary').textContent =
        checked.length === 0 ? 'All' : checked.length === 1 ? checked[0] : `${checked.length} selected`;
      ms.classList.toggle('active', checked.length > 0);
    });
  }

  // Charts and the demux table (both further down the page).
  function applyFilterElsewhere(f) {
    document.querySelectorAll('#demux-table tbody.run-group').forEach(tb => {
      tb.style.display = matchesFilter(f, tb) ? '' : 'none';
    });
    document.querySelectorAll('.plot-container iframe').forEach(fr => {
      try { fr.contentWindow.setFilter && fr.contentWindow.setFilter(f); } catch (e) {}
    });
  }

  function applyFilters() {
    const f = currentTableFilter();
    updateFilterLabels();
    const rows = document.querySelectorAll('#cell-table ' + CELL_ROWS);
    let shown = 0;
    rows.forEach(r => {
      const match = matchesFilter(f, r);
      r.style.display = match ? '' : 'none';
      if (match) shown++;
    });
    document.getElementById('row-count').textContent = `Showing ${shown} of ${rows.length} acquisitions`;
    updateSums();
    applyDataBars();
    applyFilterElsewhere(f);
  }

  const COLOR_ATTR = { sample: 'colorSample', run: 'colorRun', cell: 'colorCell' };
  const NAME_ATTR  = { sample: 'sample',       run: 'run',       cell: 'cell' };

  const TAG_COL = { sample: 'sample_name', run: 'run_name', cell: 'cell_id' };

  // Stripe plus a tinted tag on the column that names the colour group.
  function applyRowColors(mode) {
    const ci = HEADER_THS.findIndex(th => th.dataset.key === TAG_COL[mode]);
    document.querySelectorAll('#cell-table ' + CELL_ROWS).forEach(r => {
      const color = r.dataset[COLOR_ATTR[mode]];
      const stripe = r.querySelector('td.color-stripe');
      if (stripe) stripe.style.boxShadow = `inset 4px 0 0 ${color},inset -1px 0 0 #e2e8f0`;
      r.querySelectorAll('.tag').forEach(t => {
        const on = t.parentElement === r.children[ci];
        t.classList.toggle('on', on);
        t.style.background = on ? `color-mix(in srgb, ${color} 16%, transparent)` : '';
        t.style.color = on ? `color-mix(in srgb, ${color} 75%, #000)` : '';
      });
    });
  }

  function buildTableLegend(mode) {
    const seen = new Map();
    document.querySelectorAll('#cell-table ' + CELL_ROWS).forEach(r => {
      const name = r.dataset[NAME_ATTR[mode]] || '(none)';
      if (!seen.has(name)) seen.set(name, r.dataset[COLOR_ATTR[mode]]);
    });
    const legend = document.getElementById('table-legend');
    legend.innerHTML = '';
    seen.forEach((color, name) => {
      const item = document.createElement('div');
      item.style.display = 'flex';
      item.style.alignItems = 'center';
      item.style.gap = '5px';
      item.innerHTML = `<span style="width:10px;height:10px;border-radius:3px;` +
        `background:${color};display:inline-block;"></span>${name}`;
      legend.appendChild(item);
    });
  }

  document.getElementById('color-by-select').addEventListener('change', e => {
    applyRowColors(e.target.value);
    buildTableLegend(e.target.value);
  });
  applyRowColors('sample');
  buildTableLegend('sample');

  // ---------------------------------------------------------------------------
  // DATA BARS: thin inset pill under zero-based amounts (data-bar columns).
  // Largest visible value per column = 100% fill.
  // ---------------------------------------------------------------------------
  const BAR_COL_INDICES = HEADER_THS.map((th, i) => th.dataset.bar === 'true' ? i : -1).filter(i => i >= 0);

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
    const rows = visibleCellRows(document.getElementById('cell-table'));
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

  FILTER_IDS.forEach(id => document.getElementById(id).addEventListener('change', applyFilters));
  document.getElementById('filter-clear').addEventListener('click', () => {
    document.querySelectorAll('.ms input:checked').forEach(i => { i.checked = false; });
    applyFilters();
  });
  // One dropdown open at a time; a click elsewhere closes it.
  document.querySelectorAll('.ms').forEach(ms => ms.addEventListener('toggle', () => {
    if (ms.open) document.querySelectorAll('.ms[open]').forEach(o => { if (o !== ms) o.open = false; });
  }));
  document.addEventListener('click', e => {
    document.querySelectorAll('.ms[open]').forEach(o => { if (!o.contains(e.target)) o.open = false; });
  });
  applyFilters();
  // The demux table and later charts are parsed after this script.
  document.addEventListener('DOMContentLoaded', () => applyFilterElsewhere(currentTableFilter()));
</script>

<h3>Metric by Acquisition</h3>
<div class="plot-container">
  <iframe srcdoc="%%BAR_IFRAME%%" height="600" frameborder="0" scrolling="no"></iframe>
</div>

%%REUSE_SECTION%%
%%DEMUX_SECTION%%
<h3>Metric Correlation</h3>
<div class="plot-container">
  <iframe srcdoc="%%SCATTER_IFRAME%%" height="600" frameborder="0" scrolling="no"></iframe>
</div>

%%FOOTER%%
%%PAGE_TIP%%
</body>
</html>
"""


# ---------------------------------------------------------------------------
# IFRAMES: shared styles and JS (on top of report_common.CHART_HEAD)
# ---------------------------------------------------------------------------

_IFRAME_EXTRA = """\
  <style>
    h2{font-size:1.1rem;}
    .controls{display:flex;justify-content:center;align-items:center;gap:10px;
               margin-bottom:14px;flex-wrap:wrap;}
    .controls label{font-size:12.5px;color:#555;font-weight:600;}
    .controls select{
      padding:5px 10px;font-size:12.5px;border:1px solid #cbd5e0;border-radius:6px;
      background:#f8fafc;color:#2d3748;cursor:pointer;max-width:100%;
    }
    .bar{opacity:.9;}
    .bar:hover{opacity:1;}
    .dot{opacity:.85;stroke:#fff;stroke-width:1px;}
    .dot:hover{opacity:1;}
  </style>
  <script>
    window.MAX_CHART_W = 1020;

    // Group legend (by sample / run / SMRT Cell): right of the plot when there
    // is room, else in rows below. Click hides a group; build() redraws.
    function groupLegendLayout(totalW, n) {
      if (totalW >= NARROW) return {right: 160, below: false, cols: 1, height: 0};
      const cols = Math.max(1, Math.floor((totalW - 30) / 180));
      return {right: 16, below: true, cols, height: Math.ceil(n / cols) * 22 + 8};
    }

    // The report's Acquisition Summary filters: {run, cell, sample} lists,
    // empty = all. Read from the parent on start, then pushed via setFilter().
    let tableFilter = {run: [], cell: [], sample: []};
    try { if (parent.currentTableFilter) tableFilter = parent.currentTableFilter(); } catch (e) {}
    // Before the chart's own script has run, build() will pick the filter up.
    window.setFilter = f => { tableFilter = f; if (document.readyState !== 'loading') build(); };
    // Works for chart records and demux groups (both have run, cellId, sample).
    const inFilter = r => (!tableFilter.run.length || tableFilter.run.includes(r.run)) &&
                          (!tableFilter.cell.length || tableFilter.cell.includes(r.cellId)) &&
                          (!tableFilter.sample.length || tableFilter.sample.includes(r.sample));
    // Legend entries for the groups that pass the filters.
    function legendGroups() {
      const present = new Set(records.filter(inFilter).map(groupNameOf));
      return (groups[colorMode] || []).filter(g => present.has(g.name));
    }
    function noDataNote() {
      d3.select('#chart').append('p').style('text-align', 'center').style('color', '#718096')
        .style('padding', '40px 0').text('No acquisitions match the filters.');
      fitFrame();
    }

    function groupNameOf(r) {
      if (colorMode === 'run') return r.run || 'Unknown';
      if (colorMode === 'cell') return r.cellId || 'Unknown';
      return r.sample || 'Unknown';
    }

    function drawGroupLegend(svg, W, H, legend, axisSpace, shape) {
      const g = svg.append('g').attr('transform',
        legend.below ? `translate(${-40},${H + axisSpace})` : `translate(${W + 14},0)`);
      legendGroups().forEach((grp, i) => {
        const hidden = hiddenGroups.has(grp.name);
        const pos = legend.below
          ? [(i % legend.cols) * 180, Math.floor(i / legend.cols) * 22 + 6]
          : [0, i * 20];
        const row = g.append('g').attr('class', 'legend-item')
            .attr('transform', `translate(${pos[0]},${pos[1]})`)
            .on('click', () => {
              if (hidden) hiddenGroups.delete(grp.name); else hiddenGroups.add(grp.name);
              build();
            });
        if (shape === 'circle') {
          row.append('circle').attr('cx', 5.5).attr('cy', 5.5).attr('r', 5.5)
             .attr('fill', hidden ? '#ccc' : grp.color);
        } else {
          row.append('rect').attr('width', 11).attr('height', 11).attr('rx', 2)
             .attr('fill', hidden ? '#ccc' : grp.color);
        }
        row.append('text').attr('class', 'legend-text').attr('x', 16).attr('y', 9)
           .attr('fill', hidden ? '#999' : '#333')
           .style('text-decoration', hidden ? 'line-through' : 'none')
           .text(grp.name);
      });
    }

    // Bar labels for acquisitions: the well, then "SMRT Cell · use N"; under
    // each run's bars a bracket with the run name. Falls back to one rotated
    // line per bar when bars are too narrow for two flat lines.
    function drawAcquisitionAxis(svg, data, x, H) {
      const g = svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`);
      const mid = d => x(d.label) + x.bandwidth() / 2;
      const use = d => d.metrics.acquisition_no ? ` · use ${d.metrics.acquisition_no}` : '';
      if (x.step() >= 84) {
        data.forEach(d => {
          const t = g.append('text').attr('text-anchor', 'middle').attr('y', 16);
          t.append('tspan').attr('x', mid(d)).style('font-weight', '600').text(d.well || d.label);
          t.append('tspan').attr('x', mid(d)).attr('dy', '1.25em').style('font-size', '10.5px')
           .style('fill', '#718096').text(`${d.cellId}${use(d)}`);
        });
      } else {
        data.forEach(d => {
          g.append('text').attr('transform', `translate(${mid(d)},10) rotate(-60)`)
           .attr('text-anchor', 'end').style('font-size', '10.5px')
           .text(`${d.well} · ${d.cellId}${use(d)}`);
        });
      }
      // Runs: consecutive bars of the same run share a bracket.
      const runs = [];
      data.forEach(d => {
        const last = runs[runs.length - 1];
        if (last && last.run === d.run) last.items.push(d);
        else runs.push({run: d.run, items: [d]});
      });
      // Brackets go just below the bar labels, however far those reach.
      const by = g.node().getBBox().y + g.node().getBBox().height + 10;
      runs.forEach(({run, items}) => {
        const x0 = x(items[0].label), x1 = x(items[items.length - 1].label) + x.bandwidth();
        g.append('path').attr('d', `M${x0},${by - 5}V${by}H${x1}V${by - 5}`)
         .attr('fill', 'none').attr('stroke', '#a0aec0');
        g.append('text').attr('x', (x0 + x1) / 2).attr('y', by + 15).attr('text-anchor', 'middle')
         .style('fill', '#4a5568').text(run || 'Unknown run');
      });
      return by + 22;   // height used below the plot
    }

    function fillSelect(sel, options) {
      sel.selectAll('option').data(options).join('option')
         .attr('value', d => d.key).text(d => d.label);
    }

    // Tooltip for one acquisition: highlighted value rows first, then identity.
    function cellTip(d, valueRows) {
      return tipHtml(escTip(d.label), [
        ...valueRows,
        [null, 'Sample', escTip(d.sample)],
        [null, 'Run', escTip(d.run)],
        [null, 'SMRT Cell', `${escTip(d.cellId)}${d.metrics.acquisition_no ? ` (use ${d.metrics.acquisition_no})` : ''}`],
      ]);
    }
  </script>"""


# ---------------------------------------------------------------------------
# IFRAME 1: BAR CHART -- dropdown selects which metric to plot per acquisition
# ---------------------------------------------------------------------------

BAR_CHART_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>Metric by Acquisition</title>
%%CHART_HEAD%%
%%IFRAME_EXTRA%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
%%EXPORT_BUTTONS%%
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
    const payload = %%PLOT_DATA%%;
    const records = payload.records;
    const metrics = payload.metrics;
    const groups = payload.groups;
    const colorModes = payload.colorModes;
    let colorMode = colorModes.length ? colorModes[0].key : 'sample';
    let hiddenGroups = new Set();

    const select = d3.select('#metric-select');
    fillSelect(select, metrics);
    const colorSelect = d3.select('#color-by-select');
    fillSelect(colorSelect, colorModes);
    colorSelect.property('value', colorMode);

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const metricKey = select.property('value');
      const meta = metrics.find(m => m.key === metricKey) || {};
      document.getElementById('chart-title').textContent = meta.label || metricKey;

      const totalW = chartWidth();
      const legend = groupLegendLayout(totalW, legendGroups().length);
      const margin = {top: 20, right: legend.right, bottom: 84, left: 70};
      const W = totalW - margin.left - margin.right;
      const H = 360;
      const svg = appendSvg(totalW, margin, H, legend, 92);

      const data = records.filter(r => r.metrics[metricKey] !== null &&
                                       r.metrics[metricKey] !== undefined &&
                                       inFilter(r) && !hiddenGroups.has(groupNameOf(r)));
      if (!data.length) { d3.select('#chart').selectAll('*').remove(); return noDataNote(); }

      const x = d3.scaleBand().domain(data.map(d => d.label)).range([0, W]).padding(0.3);
      const y = d3.scaleLinear()
          .domain([0, (d3.max(data, d => d.metrics[metricKey]) || 1) * 1.08]).range([H, 0]);

      svg.append('g').attr('class', 'grid').call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
      svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x).tickFormat(''));
      svg.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickFormat(d3.format(',.2~f')));
      const axisH = drawAcquisitionAxis(svg, data, x, H);
      axisLabel(svg, meta.label || metricKey, {transform: 'rotate(-90)', x: -H / 2, y: -55});

      const decimals = meta.decimals ?? 2;
      svg.append('g').selectAll('rect').data(data, d => d.label)
        .join('rect')
        .attr('class', 'bar')
        .attr('x', d => x(d.label))
        .attr('width', x.bandwidth())
        .attr('fill', d => d.colors[colorMode])
        .attr('y', d => y(d.metrics[metricKey]))
        .attr('height', d => H - y(d.metrics[metricKey]))
        .on('mousemove', (event, d) => showTooltip(event, cellTip(d, [
          [d.colors[colorMode], meta.label, d3.format(',.' + decimals + 'f')(d.metrics[metricKey])]])))
        .on('mouseout', hideTooltip);

      drawGroupLegend(svg, W, H, legend, axisH + 4, 'rect');
      // Size the chart to the labels actually drawn (rotated ones need more room).
      d3.select(svg.node().ownerSVGElement)
        .attr('height', margin.top + H + axisH + (legend.below ? legend.height + 4 : 0));
      fitFrame();
    }

    select.on('change', build);
    colorSelect.on('change', function () {
      colorMode = this.value;
      hiddenGroups = new Set();
      build();
    });
    if (metrics.length) {
      select.property('value', metrics[0].key);
      build();
      onResize(build);
    }
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
  <title>Metric Correlation</title>
%%CHART_HEAD%%
%%IFRAME_EXTRA%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
%%EXPORT_BUTTONS%%
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
    const payload = %%PLOT_DATA%%;
    const records = payload.records;
    const metrics = payload.metrics;
    const groups = payload.groups;
    const colorModes = payload.colorModes;
    let colorMode = colorModes.length ? colorModes[0].key : 'sample';
    let hiddenGroups = new Set();

    const xSelect = d3.select('#x-select');
    const ySelect = d3.select('#y-select');
    fillSelect(xSelect, metrics);
    fillSelect(ySelect, metrics);
    const colorSelect = d3.select('#color-by-select');
    fillSelect(colorSelect, colorModes);
    colorSelect.property('value', colorMode);

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const xKey = xSelect.property('value');
      const yKey = ySelect.property('value');
      const xMeta = metrics.find(m => m.key === xKey) || {};
      const yMeta = metrics.find(m => m.key === yKey) || {};
      document.getElementById('chart-title').textContent =
        `${yMeta.label || yKey} vs ${xMeta.label || xKey}`;

      const totalW = chartWidth();
      const legend = groupLegendLayout(totalW, legendGroups().length);
      const margin = {top: 20, right: legend.right, bottom: 60, left: 70};
      const W = totalW - margin.left - margin.right;
      const H = 460 - margin.top - margin.bottom;
      const svg = appendSvg(totalW, margin, H, legend, 60);

      const data = records.filter(r =>
        r.metrics[xKey] !== null && r.metrics[xKey] !== undefined &&
        r.metrics[yKey] !== null && r.metrics[yKey] !== undefined &&
        inFilter(r) && !hiddenGroups.has(groupNameOf(r)));
      if (!data.length) { d3.select('#chart').selectAll('*').remove(); return noDataNote(); }

      const xExtent = d3.extent(data, d => d.metrics[xKey]);
      const yExtent = d3.extent(data, d => d.metrics[yKey]);
      const xPad = (xExtent[1] - xExtent[0]) * 0.1 || 1;
      const yPad = (yExtent[1] - yExtent[0]) * 0.1 || 1;
      const x = d3.scaleLinear().domain([xExtent[0] - xPad, xExtent[1] + xPad]).range([0, W]);
      const y = d3.scaleLinear().domain([yExtent[0] - yPad, yExtent[1] + yPad]).range([H, 0]);
      const ticks = Math.max(3, Math.min(8, Math.floor(W / 90)));

      svg.append('g').attr('class', 'grid').call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
      svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x).ticks(ticks).tickFormat(d3.format(',.2~f')));
      svg.append('g').attr('class', 'axis')
         .call(d3.axisLeft(y).ticks(8).tickFormat(d3.format(',.2~f')));
      axisLabel(svg, xMeta.label || xKey, {x: W / 2, y: H + 45});
      axisLabel(svg, yMeta.label || yKey, {transform: 'rotate(-90)', x: -H / 2, y: -55});

      const xd = xMeta.decimals ?? 2, yd = yMeta.decimals ?? 2;
      svg.append('g').selectAll('circle').data(data, d => d.label)
        .join('circle')
        .attr('class', 'dot')
        .attr('r', 6)
        .attr('fill', d => d.colors[colorMode])
        .attr('cx', d => x(d.metrics[xKey]))
        .attr('cy', d => y(d.metrics[yKey]))
        .on('mousemove', (event, d) => showTooltip(event, cellTip(d, [
          [d.colors[colorMode], yMeta.label, d3.format(',.' + yd + 'f')(d.metrics[yKey])],
          [null, xMeta.label, d3.format(',.' + xd + 'f')(d.metrics[xKey])]])))
        .on('mouseout', hideTooltip);

      drawGroupLegend(svg, W, H, legend, 60, 'circle');
      fitFrame();
    }

    xSelect.on('change', build);
    ySelect.on('change', build);
    colorSelect.on('change', function () {
      colorMode = this.value;
      hiddenGroups = new Set();
      build();
    });

    if (metrics.length) {
      xSelect.property('value', metrics.find(m => m.key === 'loading_conc_pM') ? 'loading_conc_pM' : metrics[0].key);
      ySelect.property('value', metrics.find(m => m.key === 'hifi_yield_gb') ? 'hifi_yield_gb' : (metrics[1] || metrics[0]).key);
      build();
      onResize(build);
    }
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# IFRAME 3: SAMPLE SUMMARY -- HiFi yield and reads per sample, one segment per acquisition
# ---------------------------------------------------------------------------

SAMPLE_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>HiFi Yield and Reads per Sample</title>
%%CHART_HEAD%%
%%IFRAME_EXTRA%%
  <style>
    .export-bar{position:static;justify-content:flex-end;margin-bottom:2px;}
    #chart-container{max-width:none;}   /* use the section's full width */
  </style>
  <script>window.MAX_CHART_W = 4000;</script>
</head>
<body>
  <div id="chart-container">
    <!-- The section heading names this chart; the hidden title is used for exports. -->
    <h2 id="chart-title" style="display:none">HiFi Yield and Reads per Sample</h2>
%%EXPORT_BUTTONS%%
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const payload = %%PLOT_DATA%%;
    let colorMode = 'sample';

    // Two panels sharing the sample rows: HiFi yield and HiFi (Q20+) reads.
    const fmtReads = n => n >= 1e6 ? (n / 1e6).toFixed(1) + ' M' : n.toLocaleString();
    const PANELS = [
      {key: 'hifi_yield_gb', title: 'HiFi yield (Gb)', label: 'HiFi yield', share: 'Share of sample yield',
       fmt: v => v.toFixed(1) + ' Gb', tipFmt: v => v.toFixed(2) + ' Gb'},
      {key: 'q20_reads', title: 'HiFi reads (Q20+)', label: 'HiFi reads', share: 'Share of sample reads',
       fmt: fmtReads, tipFmt: v => Math.round(v).toLocaleString()},
    ];

    // One bar per sample and panel; each acquisition (one use of a SMRT
    // Cell) is a segment, in run order, with its [x0, x1] per panel.
    function aggregate(records) {
      const samples = [];
      const bySample = new Map();
      records.forEach(r => {
        const name = r.sample || 'Unknown';
        if (!bySample.has(name)) {
          const entry = {name, color: r.colors.sample, cells: [], totals: {}};
          PANELS.forEach(p => { entry.totals[p.key] = 0; });
          bySample.set(name, entry);
          samples.push(entry);
        }
        const entry = bySample.get(name);
        const seg = {};
        PANELS.forEach(p => {
          const v = r.metrics[p.key] || 0;
          seg[p.key] = [entry.totals[p.key], entry.totals[p.key] + v];
          entry.totals[p.key] += v;
        });
        entry.cells.push({r, seg});
      });
      return samples;
    }

    function segmentTip(s, d, panel) {
      const [x0, x1] = d.seg[panel.key];
      return cellTip(d.r, [
        ...PANELS.map((p, i) => [i === 0 ? s.color : null, p.label,
                                 p.tipFmt(d.seg[p.key][1] - d.seg[p.key][0])]),
        [null, panel.share,
         `${(100 * (x1 - x0) / (s.totals[panel.key] || 1)).toFixed(0)}%`],
      ]);
    }

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const samples = aggregate(payload.records.filter(inFilter));
      if (!samples.length) return noDataNote();
      const totalW = chartWidth();
      const narrow = totalW < NARROW;   // panels stacked instead of side by side
      const labelW = 136, gap = 44, valueW = 74, titleH = 22, axisH = 30;
      const H = samples.length * 36;
      const blockH = titleH + H + axisH;
      const panelW = narrow ? totalW - labelW : (totalW - labelW - gap) / 2;
      const W = Math.max(60, panelW - valueW);
      const root = d3.select('#chart').append('svg')
          .attr('width', totalW).attr('height', narrow ? 2 * blockH + 12 : blockH);
      const y = d3.scaleBand().domain(samples.map(s => s.name)).range([0, H]).padding(0.25);

      PANELS.forEach((panel, pi) => {
        const ox = narrow ? labelW : labelW + pi * (panelW + gap);
        const oy = (narrow ? pi * (blockH + 12) : 0) + titleH;
        const g = root.append('g').attr('transform', `translate(${ox},${oy})`);
        const x = d3.scaleLinear()
            .domain([0, (d3.max(samples, s => s.totals[panel.key]) || 1) * 1.02]).range([0, W]);
        const ticks = Math.max(2, Math.floor(W / 70));

        g.append('text').attr('x', 0).attr('y', -8)
         .style('font-size', '12.5px').style('font-weight', '600').style('fill', '#333')
         .text(panel.title);
        g.append('g').attr('class', 'grid').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x).ticks(ticks).tickSize(-H).tickFormat(''));
        g.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x).ticks(ticks).tickFormat(panel.key === 'q20_reads' ? d3.format('.2~s') : null));

        // Sample names (and acquisition counts) once per row of panels.
        if (pi === 0 || narrow) {
          samples.forEach(s => {
            const n = s.cells.length;
            const smrt = new Set(s.cells.map(c => c.r.cellId)).size;
            const t = g.append('text').attr('x', -10).attr('text-anchor', 'end')
              .attr('y', y(s.name) + y.bandwidth() / 2).style('font-size', '12px').style('fill', '#333');
            t.append('tspan').attr('dy', '-0.1em').text(s.name);
            t.append('tspan').attr('x', -10).attr('dy', '1.15em')
             .style('font-size', '10.5px').style('fill', '#718096')
             .text(`${n} acq. · ${smrt} SMRT Cell${smrt === 1 ? '' : 's'}`);
          });
        }

        samples.forEach(s => {
          const row = g.append('g');
          row.selectAll('rect').data(s.cells).join('rect')
            .attr('class', 'bar')
            .attr('x', d => x(d.seg[panel.key][0]))
            .attr('y', y(s.name))
            .attr('width', d => Math.max(0, x(d.seg[panel.key][1]) - x(d.seg[panel.key][0])))
            .attr('height', y.bandwidth())
            // Alternate shades so neighbouring acquisitions stay distinguishable.
            .attr('fill', (d, i) => i % 2 ? d3.color(s.color).brighter(0.55) : s.color)
            .attr('stroke', '#fff').attr('stroke-width', 1)
            .on('mousemove', (event, d) => showTooltip(event, segmentTip(s, d, panel)))
            .on('mouseout', hideTooltip);
          row.append('text').attr('class', 'legend-text')
            .attr('x', x(s.totals[panel.key]) + 8).attr('y', y(s.name) + y.bandwidth() / 2)
            .attr('dy', '0.35em').attr('fill', '#333').style('font-weight', '600')
            .text(panel.fmt(s.totals[panel.key]));
        });
      });
      fitFrame();
    }

    build();
    onResize(build);
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# IFRAME 4: MULTI-USE CELLS -- a metric across acquisitions of the same cell
# ---------------------------------------------------------------------------

REUSE_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>Multi-use Cells</title>
%%CHART_HEAD%%
%%IFRAME_EXTRA%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
%%EXPORT_BUTTONS%%
    <div class="controls">
      <label for="metric-select">Metric</label>
      <select id="metric-select"></select>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const payload = %%PLOT_DATA%%;
    const metrics = payload.metrics;
    // Only cells acquired more than once; lines and legend are per cell.
    const counts = d3.rollup(payload.records, v => v.length, r => r.cellId);
    const records = payload.records.filter(r => counts.get(r.cellId) > 1);
    const colorMode = 'cell';
    const groups = {cell: payload.groups.cell.filter(g => counts.get(g.name) > 1)};
    let hiddenGroups = new Set();

    const select = d3.select('#metric-select');
    fillSelect(select, metrics.filter(m => m.key !== 'acquisition_no'));

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const key = select.property('value');
      const meta = metrics.find(m => m.key === key) || {};
      document.getElementById('chart-title').textContent =
        `${meta.label || key} by Acquisition (multi-use cells)`;

      const totalW = chartWidth();
      const legend = groupLegendLayout(totalW, legendGroups().length);
      const margin = {top: 20, right: legend.right, bottom: 50, left: 70};
      const W = totalW - margin.left - margin.right;
      const H = 380 - margin.top - margin.bottom;
      const svg = appendSvg(totalW, margin, H, legend, 56);

      const data = records.filter(r => r.metrics[key] != null && r.metrics.acquisition_no != null &&
                                       inFilter(r) && !hiddenGroups.has(r.cellId));
      if (!data.length) { d3.select('#chart').selectAll('*').remove(); return noDataNote(); }
      const lines = d3.groups(data, r => r.cellId)
        .map(([cell, rs]) => [cell, rs.sort((a, b) => a.metrics.acquisition_no - b.metrics.acquisition_no)]);
      const maxAcq = d3.max(records, r => r.metrics.acquisition_no) || 1;
      const x = d3.scaleLinear().domain([0.8, maxAcq + 0.2]).range([0, W]);
      const ext = d3.extent(data, r => r.metrics[key]);
      const pad = (ext[1] - ext[0]) * 0.15 || Math.abs(ext[1] || 1) * 0.1;
      const y = d3.scaleLinear().domain([ext[0] - pad, ext[1] + pad]).nice().range([H, 0]);

      svg.append('g').attr('class', 'grid').call(d3.axisLeft(y).tickSize(-W).tickFormat(''));
      svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${H})`)
         .call(d3.axisBottom(x).tickValues(d3.range(1, maxAcq + 1)).tickFormat(d3.format('d')));
      svg.append('g').attr('class', 'axis').call(d3.axisLeft(y).tickFormat(d3.format(',.2~f')));
      axisLabel(svg, 'Acquisition', {x: W / 2, y: H + 42});
      axisLabel(svg, meta.label || key, {transform: 'rotate(-90)', x: -H / 2, y: -55});

      const line = d3.line().x(r => x(r.metrics.acquisition_no)).y(r => y(r.metrics[key]));
      const decimals = meta.decimals ?? 2;
      lines.forEach(([cell, rs]) => {
        const color = rs[0].colors.cell;
        svg.append('path').datum(rs).attr('fill', 'none').attr('stroke', color)
           .attr('stroke-width', 2.5).attr('opacity', 0.85).attr('d', line);
        svg.append('g').selectAll('circle').data(rs).join('circle')
           .attr('class', 'dot').attr('r', 5.5).attr('fill', color)
           .attr('cx', r => x(r.metrics.acquisition_no)).attr('cy', r => y(r.metrics[key]))
           .on('mousemove', (event, r) => showTooltip(event, cellTip(r, [
             [color, meta.label, d3.format(',.' + decimals + 'f')(r.metrics[key])]])))
           .on('mouseout', hideTooltip);
      });

      drawGroupLegend(svg, W, H, legend, 56, 'circle');
      fitFrame();
    }

    select.on('change', build);
    select.property('value', metrics.find(m => m.key === 'hifi_yield_gb') ? 'hifi_yield_gb' : metrics[0].key);
    build();
    onResize(build);
  </script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# IFRAME 5: DEMULTIPLEXING -- a metric per barcode, grouped by acquisition
# ---------------------------------------------------------------------------

DEMUX_IFRAME = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <title>Demultiplexing</title>
%%CHART_HEAD%%
%%IFRAME_EXTRA%%
</head>
<body>
  <div id="chart-container">
    <h2 id="chart-title"></h2>
%%EXPORT_BUTTONS%%
    <div class="controls">
      <label for="metric-select">Metric</label>
      <select id="metric-select"></select>
    </div>
    <div id="chart"></div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <script>
    const payload = %%PLOT_DATA%%;
    const groups = payload.groups;
    const metrics = payload.metrics;
    const select = d3.select('#metric-select');
    fillSelect(select, metrics);
    // Set from the table's "Color by" (see applyDemuxColors in the report).
    let colorMode = 'sample';
    window.setColorMode = mode => { if (mode !== colorMode) { colorMode = mode; build(); } };

    function build() {
      d3.select('#chart').selectAll('*').remove();
      const key = select.property('value');
      const meta = metrics.find(m => m.key === key) || {};
      document.getElementById('chart-title').textContent = `${meta.label || key} per Barcode`;
      const fmt = d3.format(',.' + (meta.decimals ?? 2) + 'f');
      const shown = groups.filter(inFilter);
      if (!shown.length) return noDataNote();

      const totalW = chartWidth();
      const labelW = Math.min(170, totalW * 0.3), valueW = 70, rowH = 24, headH = 34, axisH = 26;
      const W = Math.max(60, totalW - labelW - valueW);
      const H = shown.reduce((h, g) => h + headH + g.rows.length * rowH + axisH, 0);
      const svg = d3.select('#chart').append('svg').attr('width', totalW).attr('height', H)
        .append('g').attr('transform', `translate(${labelW},0)`);
      const ticks = Math.max(2, Math.floor(W / 80));

      let oy = 0;
      shown.forEach(g => {
        const head = svg.append('text').attr('x', -labelW).attr('y', oy + headH - 10)
          .style('font-size', '12.5px').style('fill', '#333');
        head.append('tspan').style('font-weight', '600').text(g.label);
        if (g.sample) head.append('tspan').style('fill', '#718096').text(`  ${g.sample} · ${g.cellId}`);
        oy += headH;

        // Each acquisition has its own scale, to show how evenly it split.
        const gH = g.rows.length * rowH;
        const max = d3.max(g.rows, r => r.metrics[key]) || 1;
        const x = d3.scaleLinear().domain([0, max * 1.04]).range([0, W]);
        svg.append('g').attr('class', 'grid').attr('transform', `translate(0,${oy + gH})`)
           .call(d3.axisBottom(x).ticks(ticks).tickSize(-gH).tickFormat(''));
        svg.append('g').attr('class', 'axis').attr('transform', `translate(0,${oy + gH})`)
           .call(d3.axisBottom(x).ticks(ticks).tickFormat(d3.format('.3~s')));

        // Even split: what each barcode would get if reads were spread equally.
        const bc = g.rows.filter(r => r.barcoded && r.metrics[key] != null);
        if (meta.even && bc.length > 1) {
          const even = d3.sum(bc, r => r.metrics[key]) / bc.length;
          svg.append('line').attr('x1', x(even)).attr('x2', x(even))
             .attr('y1', oy - 2).attr('y2', oy + bc.length * rowH + 2)
             .attr('stroke', '#475569').attr('stroke-dasharray', '4,3')
             .append('title').text(`Even split: ${fmt(even)}`);
        }

        g.rows.forEach((r, i) => {
          const y = oy + i * rowH, v = r.metrics[key];
          svg.append('text').attr('x', -8).attr('y', y + rowH / 2).attr('dy', '0.35em')
             .attr('text-anchor', 'end').style('font-size', '11.5px')
             .style('fill', r.barcoded ? '#333' : '#94a3b8').text(r.name);
          if (v == null) return;
          svg.append('rect').attr('class', 'bar').attr('x', 0).attr('y', y + 4)
             .attr('width', x(v)).attr('height', rowH - 8).attr('rx', 2)
             .attr('fill', r.barcoded ? r.colors[colorMode] : '#cbd5e0')
             .on('mousemove', event => showTooltip(event, tipHtml(escTip(r.name), [
               [r.barcoded ? r.colors[colorMode] : '#cbd5e0', meta.label, fmt(v)],
               [null, 'Barcode', escTip(r.barcode)],
               [null, 'Acquisition', escTip(g.label)],
               [null, 'Sample', escTip(g.sample || '-')],
             ])))
             .on('mouseout', hideTooltip);
          svg.append('text').attr('class', 'legend-text').attr('x', x(v) + 6)
             .attr('y', y + rowH / 2).attr('dy', '0.35em').attr('fill', '#4a5568')
             .text(key === 'reads' && v >= 1e6 ? (v / 1e6).toFixed(2) + ' M' : fmt(v));
        });
        oy += gH + axisH;
      });
      fitFrame();
    }

    select.on('change', build);
    select.property('value', metrics[0].key);
    build();
    onResize(build);
  </script>
</body>
</html>"""


def _render_iframe(template, payload, title):
    return render_chart_iframe(template, payload, title, iframe_extra=_IFRAME_EXTRA)


# ---------------------------------------------------------------------------
# MAIN GENERATION FUNCTION
# ---------------------------------------------------------------------------

def generate_report(files, title, out_path, demux_files=()):
    cells = load_all_cells(files)
    if not cells:
        sys.exit("ERROR: No SMRT Cell rows found in the given file(s).")

    assign_colors(cells)
    demuxes = load_all_demux(demux_files, cells)

    n_runs = len({c["run_id"] or c["run_name"] for c in cells})
    n_smrt = len({c["cell_id"] for c in cells if c["cell_id"]})
    n_samp = len({c["sample_name"] or NO_SAMPLE for c in cells})
    def plural(n, word):
        return f"{n} {word}{'' if n == 1 else 's'}"
    # A row is an acquisition: one use of a SMRT Cell (a well in a run).
    counts = (f"{plural(len(cells), 'acquisition')} of {plural(n_smrt, 'SMRT Cell')} "
              f"&middot; {plural(n_samp, 'sample')} &middot; {plural(n_runs, 'run')}")
    print(f"Loaded {plural(len(cells), 'acquisition')} of {plural(n_smrt, 'SMRT Cell')} "
          f"({plural(n_samp, 'sample')}, {plural(n_runs, 'run')}) from {plural(len(files), 'file')}.")

    payload = build_chart_payload(cells)
    bar_iframe = _render_iframe(BAR_CHART_IFRAME, payload, "Metric by Acquisition")
    scatter_iframe = _render_iframe(SCATTER_IFRAME, payload, "Metric Correlation")
    sample_iframe = _render_iframe(SAMPLE_IFRAME, payload, "HiFi Yield and Reads per Sample")
    n_samples = len({rec["sample_name"] or NO_SAMPLE for rec in cells})

    # The multi-use chart only appears when some cell was acquired more than once.
    cell_uses = {}
    for rec in cells:
        cell_uses[rec["cell_id"]] = cell_uses.get(rec["cell_id"], 0) + 1
    reuse_section = ""
    if any(n > 1 for cid, n in cell_uses.items() if cid):
        reuse_iframe = _render_iframe(REUSE_IFRAME, payload, "Multi-use Cells")
        reuse_section = (
            "<h3>Multi-use Cells</h3>\n<div class=\"plot-container\">\n"
            f"  <iframe srcdoc=\"{reuse_iframe}\" height=\"520\" frameborder=\"0\" "
            "scrolling=\"no\"></iframe>\n</div>\n")


    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M")
    # d3 goes last: it is large and must not be searched for later tokens.
    report = fill(
        HTML_TEMPLATE,
        title=title,
        table_csv_js=TABLE_CSV_JS,
        counts=counts,
        n_runs=n_runs,
        table_header=render_table_header(),
        table_bodies=render_table_body(cells),
        read_length_guide_css=read_length_guide_css(_max_read_bp(cells)),
        table_footer=render_table_footer(),
        table_legend=render_table_legend(),
        sample_iframe=sample_iframe,
        sample_iframe_h=170 + 36 * n_samples,
        reuse_section=reuse_section,
        demux_section=render_demux_section(demuxes),
        filter_run_options=render_filter_options(cells, "run_name"),
        filter_cell_options=render_filter_options(cells, "cell_id"),
        filter_sample_options=render_filter_options(cells, "sample_name"),
        bar_iframe=bar_iframe,
        scatter_iframe=scatter_iframe,
        generated_on=generated_on,
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
        description="Generate an interactive PacBio HiFi run QC HTML report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="run-qc-export.csv file(s), and optionally barcodes.summary.csv "
                             "file(s) (recognised by their columns). Auto-discovered if omitted.")
    parser.add_argument("--demux", nargs="+", default=[], metavar="FILE",
                        help="SMRT Link *.barcodes.summary.csv file(s) for the Demultiplexing section")
    parser.add_argument("--title", default="HiFi Sequencing Report",
                        help="Report title (default: 'HiFi Sequencing Report')")
    parser.add_argument("--out", default="hifi-report.html",
                        help="Output HTML file (default: hifi-report.html)")
    args = parser.parse_args()

    given = args.files or (glob.glob("*qc-export*.csv") + glob.glob("*barcodes.summary*.csv"))
    demux_files = list(args.demux) + [f for f in given if is_demux_csv(f)]
    files = sorted(f for f in given if f not in demux_files)
    if not files:
        sys.exit("ERROR: No input files found. Pass file(s) explicitly or place "
                 "*qc-export*.csv files in the current directory.")

    generate_report(files=files, title=args.title, out_path=args.out,
                    demux_files=sorted(demux_files))


if __name__ == "__main__":
    main()
