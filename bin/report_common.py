"""
report_common.py
================
Pieces shared by generate-ont-report.py and generate-hifi-report.py: the d3
embedding, the chart-iframe <head> (styles plus JS helpers for responsive
layout, legends, tooltips and PNG/SVG export), table CSV export, the footer
and template filling. Must sit next to the two scripts.
"""

import html as html_lib
import json
from pathlib import Path

REPO_URL = "https://github.com/angelovangel/ont-run-report"

D3_PATH = Path(__file__).resolve().parent / "d3.v7.9.0.min.js"
D3_CDN = "https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"


def load_d3():
    """d3 source to embed in a report, or "" (charts then use the CDN)."""
    if D3_PATH.is_file():
        return D3_PATH.read_text(encoding="utf-8")
    print(f"WARNING: {D3_PATH} not found; charts will load d3 from {D3_CDN}")
    return ""


# Goes in the report's <head>; CHART_HEAD runs it inside each chart iframe.
D3_SOURCE_TAG = """\
<!-- d3 source, run inside each chart iframe (see CHART_HEAD) -->
<script type="text/plain" id="d3-src">%%D3_JS%%</script>"""


# Shared <head> content for every chart iframe: d3, base styles and helpers
# for responsive layout, legends, tooltips and image export.
CHART_HEAD = """\
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <script>
    // The d3 source is embedded once in the parent report. srcdoc iframes share
    // its origin, so run that source here — each chart needs its own d3, bound
    // to this document. Fall back to the CDN when opened on its own.
    (function () {
      let src;
      try { src = parent.document.getElementById('d3-src').textContent; } catch (e) {}
      if (src) {
        const el = document.createElement('script');
        el.textContent = src;
        document.head.appendChild(el);
      } else {
        document.write('<script src="%%D3_CDN%%"><\\/script>');
      }
    })();
  </script>
  <style>
    body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
         background:#fff;margin:0;padding:16px 20px 20px;color:#333;}
    #chart-container{max-width:960px;margin:0 auto;position:relative;}
    .line{fill:none;stroke-width:2.5px;opacity:0.85;}
    .tooltip{position:absolute;left:0;top:0;visibility:hidden;}
%%TIP_CSS%%
    .grid line{stroke:#ddd;stroke-dasharray:4,4;}
    .grid path{stroke-width:0;}
    .axis text{font-size:12px;fill:#555;}
    .axis path,.axis line{stroke:#bbb;}
    .axis-label{font-size:12px;fill:#555;}
    .hover-line{stroke:#888;stroke-dasharray:4,4;stroke-width:1px;pointer-events:none;opacity:0;}
    h2{text-align:center;color:#222;margin:4px 0 8px;padding:0 96px;font-size:1.25rem;}
    .legend-item{cursor:pointer;}
    .legend-text{font-size:12px;fill:#333;user-select:none;}
    .toggle-bar{display:flex;flex-wrap:wrap;justify-content:center;margin-bottom:10px;}
    .toggle-btn{padding:5px 16px;font-size:12px;cursor:pointer;border:1px solid #bbb;
      background:#f4f4f4;color:#444;transition:all .2s;user-select:none;}
    .toggle-btn:first-child{border-radius:6px 0 0 6px;}
    .toggle-btn:last-child{border-radius:0 6px 6px 0;}
    .toggle-btn.active{background:#2563eb;color:#fff;border-color:#2563eb;font-weight:600;}
    .export-bar{position:absolute;top:0;right:0;display:flex;gap:4px;}
    .export-btn{padding:3px 9px;font-size:11px;border:1px solid #ccc;background:#fff;color:#555;
      border-radius:5px;cursor:pointer;}
    .export-btn:hover{background:#f0f0f0;}
    @media (max-width:640px){
      body{padding:12px 8px;}
      h2{padding:0;font-size:1.05rem;}
      .export-bar{position:static;justify-content:center;margin-bottom:6px;}
      .toggle-btn{padding:5px 10px;}
    }
  </style>
  <script>
    const SVG_NS = 'http://www.w3.org/2000/svg';
    const NARROW = 640;

    // Inverse of _pack_series: [[x, y], ...] -> [{x, y}, ...]
    function unpack(samples, keys) {
      samples.forEach(s => keys.forEach(k => {
        if (s[k]) s[k] = s[k].map(([x, y]) => ({x, y}));
      }));
    }

    function chartWidth() {
      const el = document.getElementById('chart-container');
      return Math.max(300, Math.min(window.MAX_CHART_W || 960, el.clientWidth));
    }

    // Grow/shrink the hosting iframe to fit the chart.
    function fitFrame() {
      try {
        if (window.frameElement) {
          window.frameElement.style.height =
            Math.ceil(document.body.getBoundingClientRect().height) + 'px';
        }
      } catch (e) {}
    }

    function onResize(build) {
      let lastW = chartWidth(), timer;
      window.addEventListener('resize', () => {
        clearTimeout(timer);
        timer = setTimeout(() => {
          const w = chartWidth();
          if (w !== lastW) { lastW = w; build(); }
        }, 150);
      });
    }

    // Sample legend: right of the plot when there is room, else in rows below.
    function legendLayout(totalW, n) {
      if (totalW >= NARROW) return {right: 170, below: false, cols: 1, height: 0};
      const cols = Math.max(1, Math.floor((totalW - 30) / 140));
      return {right: 16, below: true, cols, height: Math.ceil(n / cols) * 22 + 8};
    }

    // Root <svg> sized for the plot plus (when below) axis label and legend.
    function appendSvg(totalW, margin, H, legend, axisSpace) {
      const height = margin.top + H + (legend.below ? axisSpace + legend.height : margin.bottom);
      return d3.select('#chart').append('svg').attr('width', totalW).attr('height', height)
        .append('g').attr('transform', `translate(${margin.left},${margin.top})`);
    }

    // Hover highlights a sample, click hides/shows it; hidden state lives on
    // the sample object so it survives rebuilds.
    function drawLegend(svg, samples, W, H, legend, axisSpace) {
      const g = svg.append('g').attr('transform',
        legend.below ? `translate(${-40},${H + axisSpace})` : `translate(${W + 14},0)`);
      samples.forEach((s, i) => {
        const pos = legend.below
          ? [(i % legend.cols) * 140, Math.floor(i / legend.cols) * 22 + 6]
          : [0, i * 24];
        const row = g.append('g').attr('class', 'legend-item')
            .attr('transform', `translate(${pos[0]},${pos[1]})`)
            .on('mouseover', () => {
              svg.selectAll('.sample-group').style('opacity', 0.12);
              svg.select(`.sample-${i}`).style('opacity', 1).raise();
            })
            .on('mouseout', () => svg.selectAll('.sample-group').style('opacity', 1))
            .on('click', () => { s.hidden = !s.hidden; apply(); });
        const dot = row.append('circle').attr('r', 5);
        const text = row.append('text').attr('class', 'legend-text').attr('x', 10).attr('y', 4).text(s.name);
        function apply() {
          svg.select(`.sample-${i}`).style('display', s.hidden ? 'none' : null);
          dot.style('fill', s.hidden ? '#ccc' : s.color);
          text.style('text-decoration', s.hidden ? 'line-through' : 'none')
              .style('fill', s.hidden ? '#999' : '#333');
        }
        apply();
      });
    }

    function hoursAxis(x, W) {
      return d3.axisBottom(x).ticks(Math.max(4, Math.floor(W / 60))).tickFormat(d => d + ' h');
    }

    function axisLabel(svg, text, attrs) {
      const t = svg.append('text').attr('class', 'axis-label').attr('text-anchor', 'middle').text(text);
      Object.entries(attrs).forEach(([k, v]) => t.attr(k, v));
      return t;
    }

%%TIP_JS%%

    // Tooltip positioned relative to the chart container, flipped to the left
    // of the pointer when it would overflow. Shown instantly (no fade).
    function showTooltip(event, html) {
      const tip = document.getElementById('tooltip');
      const box = document.getElementById('chart-container').getBoundingClientRect();
      tip.classList.add('tt-box');
      if (tip.innerHTML !== html) tip.innerHTML = html;
      tip.style.visibility = 'visible';
      const px = event.clientX - box.left, py = event.clientY - box.top;
      const left = px + 16 + tip.offsetWidth > box.width ? px - tip.offsetWidth - 16 : px + 16;
      const maxTop = window.innerHeight - box.top - tip.offsetHeight - 8;
      tip.style.left = Math.max(0, left) + 'px';
      tip.style.top = Math.max(0, Math.min(py - 12, maxTop)) + 'px';
    }

    function hideTooltip() {
      document.getElementById('tooltip').style.visibility = 'hidden';
    }

    function saveBlob(blob, name) {
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = name;
      document.body.appendChild(a);
      a.click();
      setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }

    // Export the chart as SVG or 2x PNG. Computed styles are inlined so the
    // file renders without this page's CSS; the chart title is added on top.
    // decorate(svg, width, height) may append extras and returns added height.
    function exportChart(format, decorate) {
      const src = document.querySelector('#chart svg');
      const clone = src.cloneNode(true);
      const props = ['fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-dasharray',
                     'stroke-opacity', 'opacity', 'font-size', 'font-family', 'font-weight',
                     'text-anchor', 'text-decoration', 'display'];
      const srcEls = src.querySelectorAll('*'), dstEls = clone.querySelectorAll('*');
      srcEls.forEach((el, i) => {
        const cs = getComputedStyle(el);
        dstEls[i].setAttribute('style', props.map(p => `${p}:${cs.getPropertyValue(p)}`).join(';'));
      });
      clone.querySelectorAll('.hover-line').forEach(el => el.remove());

      const title = document.getElementById('chart-title').textContent;
      const w = +src.getAttribute('width'), titleH = 34;
      let h = +src.getAttribute('height') + titleH;
      const body = document.createElementNS(SVG_NS, 'g');
      body.setAttribute('transform', `translate(0,${titleH})`);
      while (clone.firstChild) body.appendChild(clone.firstChild);
      const bg = document.createElementNS(SVG_NS, 'rect');
      bg.setAttribute('fill', '#fff');
      bg.setAttribute('width', w);
      const heading = document.createElementNS(SVG_NS, 'text');
      heading.setAttribute('x', w / 2);
      heading.setAttribute('y', 22);
      heading.setAttribute('style', 'text-anchor:middle;font-size:16px;font-weight:600;fill:#222;'
        + 'font-family:' + getComputedStyle(document.body).fontFamily);
      heading.textContent = title;
      clone.append(bg, heading, body);
      if (decorate) h += decorate(clone, w, h) || 0;
      bg.setAttribute('height', h);
      clone.setAttribute('xmlns', SVG_NS);
      clone.setAttribute('width', w);
      clone.setAttribute('height', h);
      clone.setAttribute('viewBox', `0 0 ${w} ${h}`);

      const xml = new XMLSerializer().serializeToString(clone);
      const name = title.replace(/[^A-Za-z0-9]+/g, '_').replace(/^_|_$/g, '') || 'chart';
      if (format === 'svg') {
        saveBlob(new Blob([xml], {type: 'image/svg+xml'}), name + '.svg');
        return;
      }
      const img = new Image();
      img.onload = () => {
        const canvas = document.createElement('canvas');
        canvas.width = w * 2;
        canvas.height = h * 2;
        const ctx = canvas.getContext('2d');
        ctx.scale(2, 2);
        ctx.drawImage(img, 0, 0);
        canvas.toBlob(blob => saveBlob(blob, name + '.png'), 'image/png');
      };
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(xml);
    }
  </script>"""


EXPORT_BUTTONS = """\
    <div class="export-bar">
      <button class="export-btn" onclick="exportChart('png'%%DECORATE%%)">PNG</button>
      <button class="export-btn" onclick="exportChart('svg'%%DECORATE%%)">SVG</button>
    </div>"""


# Tooltip look shared by chart tooltips (.tooltip.tt-box) and page hovers
# (.page-tip): dark box, bold title, rows of [swatch] label .... value.
TIP_CSS = """\
    .tt-box{background:#1e293b;color:#e2e8f0;font-size:12.5px;line-height:1.45;
            padding:8px 11px;border-radius:7px;box-shadow:0 8px 20px rgba(15,23,42,.28);
            pointer-events:none;z-index:20;white-space:nowrap;text-align:left;
            font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
    .tt-box .tt-title{font-weight:650;color:#fff;font-size:13px;}
    .tt-box .tt-text{color:#cbd5e1;white-space:normal;max-width:300px;margin-top:2px;}
    .tt-box table{border-collapse:collapse;margin:4px 0 0;width:auto;background:none;
                  box-shadow:none;font-size:12.5px;}
    .tt-box tr:hover td{background:none !important;}
    .tt-box td{padding:1px 0 !important;border:0 !important;background:none !important;
               position:static !important;font-size:12.5px;color:#e2e8f0;box-shadow:none !important;}
    .tt-box td.tt-v{text-align:right;padding-left:16px !important;font-weight:650;color:#fff;
                    font-variant-numeric:tabular-nums;}
    .tt-box .tt-sw{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:7px;
                   vertical-align:0;}"""

# tipHtml(title, rows, text): tooltip markup. rows = [[colour|null, label, value], ...];
# labels are escaped (they can be sample names), title and values are trusted.
TIP_JS = """\
    function escTip(s) {
      return String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c]);
    }
    function tipHtml(title, rows, text) {
      const body = (rows || []).map(([color, label, value]) =>
        `<tr><td>${color ? `<span class="tt-sw" style="background:${color}"></span>` : ''}`
        + `${escTip(label)}</td><td class="tt-v">${value}</td></tr>`).join('');
      return (title ? `<div class="tt-title">${title}</div>` : '')
        + (text ? `<div class="tt-text">${text}</div>` : '')
        + (body ? `<table>${body}</table>` : '');
    }"""

# Instant hovers for report pages: any element with data-tip (tooltip
# markup, see tip_html) shows it next to the pointer. Goes before </body>.
PAGE_TIP = """\
<style>
  .page-tip{position:fixed;left:0;top:0;visibility:hidden;}
""" + TIP_CSS + """
</style>
<script>
""" + TIP_JS + """
  (function () {
    const tip = document.createElement('div');
    tip.className = 'page-tip tt-box';
    document.body.appendChild(tip);
    let current = null;
    function place(e) {
      const w = tip.offsetWidth, h = tip.offsetHeight;
      const left = e.clientX + 16 + w > window.innerWidth ? e.clientX - w - 16 : e.clientX + 16;
      tip.style.left = Math.max(4, left) + 'px';
      tip.style.top = Math.max(4, Math.min(e.clientY + 14, window.innerHeight - h - 8)) + 'px';
    }
    document.addEventListener('mouseover', e => {
      const el = e.target.closest('[data-tip]');
      if (el !== current) {
        current = el;
        if (el) tip.innerHTML = el.dataset.tip;
        tip.style.visibility = el ? 'visible' : 'hidden';
      }
      if (el) place(e);
    });
    document.addEventListener('mousemove', e => { if (current) place(e); });
    document.addEventListener('scroll', () => { current = null; tip.style.visibility = 'hidden'; }, true);
  })();
</script>"""


def tip_html(title, rows=(), text=""):
    """Tooltip markup for a data-tip attribute (escape it with tip_attr).

    rows: (colour or None, label, value) tuples; a colour written "ring:#hex"
    gives an outlined round swatch. Same look as chart tooltips.
    """
    def swatch(c):
        if not c:
            return ""
        if c.startswith("ring:"):
            return (f'<span class="tt-sw" style="background:none;border-radius:50%;'
                    f'box-shadow:inset 0 0 0 2px {c[5:]}"></span>')
        return f'<span class="tt-sw" style="background:{c}"></span>'
    body = "".join(
        f'<tr><td>{swatch(c)}{html_lib.escape(str(label))}</td><td class="tt-v">{value}</td></tr>'
        for c, label, value in rows)
    return ((f'<div class="tt-title">{title}</div>' if title else "")
            + (f'<div class="tt-text">{text}</div>' if text else "")
            + (f"<table>{body}</table>" if body else ""))


def tip_attr(markup):
    """data-tip="…" attribute for tooltip markup."""
    return f'data-tip="{html_lib.escape(markup, quote=True)}"'


# downloadTableCsv(tableId, suffix): CSV of a report table's visible rows,
# subtotal rows left out. Columns come from the bottom header row (a row above
# may group them). Per column: data-csv-skip on the <th> leaves it out;
# data-csv-cols names its CSV column(s) — "A|B" splits it into A and B, taking
# "a|b" from each cell's value. A cell's value is its data-csv, else its
# data-value for number columns, else its text.
TABLE_CSV_JS = """\
  function downloadTableCsv(tableId, suffix) {
    const table = document.getElementById(tableId);
    const ths = Array.from(table.querySelector('thead tr:last-child').children);
    const esc = v => /[",\\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
    const header = [];
    ths.forEach(th => {
      if ('csvSkip' in th.dataset) return;
      if (th.dataset.csvCols) header.push(...th.dataset.csvCols.split('|'));
      else header.push(th.firstChild.textContent.trim());
    });
    const lines = [header.map(esc).join(',')];
    table.querySelectorAll('tbody tr, tfoot tr').forEach(tr => {
      if (tr.style.display === 'none' || tr.classList.contains('subtotal')) return;
      const row = [];
      ths.forEach((th, i) => {
        const td = tr.children[i];
        if ('csvSkip' in th.dataset) return;
        const value = td.dataset.csv !== undefined ? td.dataset.csv
          : th.dataset.type === 'number' ? (td.dataset.value || '') : td.textContent.trim();
        // Always one CSV field per name, so rows without split values stay aligned.
        const parts = th.dataset.csvCols ? value.split('|') : [value];
        (th.dataset.csvCols ? th.dataset.csvCols.split('|') : [0]).forEach((_, k) => row.push(parts[k] || ''));
      });
      lines.push(row.map(esc).join(','));
    });
    const blob = new Blob([lines.join('\\n') + '\\n'], {type: 'text/csv'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = (document.title.replace(/[^A-Za-z0-9]+/g, '_').replace(/^_|_$/g, '') || 'report')
      + '_' + suffix + '.csv';
    document.body.appendChild(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }"""


def footer_html(title, generated_on):
    return (
        '<div class="footer">\n'
        f'  Generated on {generated_on} &mdash; {title} &mdash;\n'
        f'  <a href="{REPO_URL}" target="_blank" rel="noopener"\n'
        f'     style="color:#64748b;text-decoration:none;">{REPO_URL.removeprefix("https://")}</a>\n'
        '</div>'
    )


def fill(template, **values):
    for key, value in values.items():
        template = template.replace(f"%%{key.upper()}%%", str(value))
    return template


def html_escape_attr(s):
    return (s.replace("&", "&amp;")
             .replace('"', "&quot;")
             .replace("<", "&lt;")
             .replace(">", "&gt;")
             .replace("'", "&#x27;"))


def script_json(data):
    """JSON safe to embed in a <script> block (labels are user-supplied)."""
    return json.dumps(data, separators=(",", ":")).replace("</", "<\\/")


def render_chart_iframe(template, data, title, decorate="", **extra):
    """Fill a chart-iframe template and escape it for a srcdoc attribute.

    Tokens: %%CHART_HEAD%%, %%EXPORT_BUTTONS%%, %%TITLE%%, %%PLOT_DATA%% plus
    any in `extra`. PLOT_DATA goes last so user-supplied labels can't inject
    other tokens.
    """
    html = fill(
        template,
        chart_head=fill(CHART_HEAD, d3_cdn=D3_CDN, tip_css=TIP_CSS, tip_js=TIP_JS),
        export_buttons=fill(EXPORT_BUTTONS, decorate=decorate),
        title=html_lib.escape(title),
        **extra,
        plot_data=script_json(data),
    )
    return html_escape_attr(html)
