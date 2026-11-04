"""One self-contained D3 report page from fixed panel contracts.

The house answer to per-figure "atomic windows": a pipeline packs its fixed
default figure set into ONE offline HTML. Each panel carries
``{id, title, caption, script, data}``; ``script`` is the same chart code the
single-figure path uses (chart-forge example or ``_shared/panels/*.js``),
mounted into the panel's own div by rewriting the canonical
``d3.select("#chart")`` token — the mount convention every chart script
follows, verified at assembly time.

Visual contract is chart-forge's own: tokens.css, micro.js, vendored d3, and
the crameri ramps are read from chart-forge (file reads only). Stdlib only —
safe to import from any skill venv.
"""
from __future__ import annotations

import html as html_lib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

FORGE = Path(__file__).resolve().parents[1] / "chart-forge"
MOUNT_TOKEN = 'd3.select("#chart")'
PROVIDER_SCRIPTS = ("templates/crameri-batlow.js", "templates/crameri-vik.js")


def example_script(name: str) -> Path:
    """Path to a chart-forge example's chart.js (the single-figure authority)."""
    path = FORGE / "examples" / name / "chart.js"
    if not path.is_file():
        raise FileNotFoundError(f"chart-forge example {name!r} has no chart.js")
    return path


def panel_script(name: str) -> Path:
    """Path to a shared report panel script (_shared/panels/<name>.js)."""
    path = Path(__file__).resolve().parent / "panels" / f"{name}.js"
    if not path.is_file():
        raise FileNotFoundError(f"report panel script {name!r} not found")
    return path


def _read(path: Path) -> str:
    return Path(path).read_text(encoding="utf-8")


def _js_safe(text: str) -> str:
    return text.replace("</script", "<\\/script")


def _d3_library() -> str:
    for cand in (FORGE / "assets" / "d3" / "vendor" / "d3.min.js",
                 FORGE.parents[1] / "reference" / "d3.min.js"):
        if cand.is_file():
            return _js_safe(_read(cand))
    raise RuntimeError("report: no vendored d3.min.js found - refusing to emit a CDN page")


def _provider(path: Path) -> str:
    """Same rewrite build.py applies to --script extras: 'export const X' -> 'const X'."""
    return _js_safe(re.sub(r"^export const ", "const ", _read(path), flags=re.M))


def json_for_html(payload) -> str:
    """Serialize JSON without allowing source text to close the script element."""
    return (json.dumps(payload, ensure_ascii=False)
            .replace("&", "\\u0026")
            .replace("<", "\\u003c")
            .replace(">", "\\u003e"))


def _section(panel: dict, *, active: bool) -> str:
    """Render one tab panel, including its skipped/error state."""
    pid = html_lib.escape(panel["id"], quote=True)
    title = html_lib.escape(panel["title"])
    status = panel.get("status", "ok")
    state_label = {"ok": "ready", "skipped": "not available", "error": "failed"}.get(
        status, status
    )
    state = (f'<span class="state state-{html_lib.escape(status, quote=True)}" '
             f'role="status">{html_lib.escape(state_label)}</span>')
    if status == "error":
        note = html_lib.escape(panel.get("note", "unknown error"))
        body = (f'<div class="state-card state-card-error"><strong>Section failed</strong>'
                f'<p>{note}</p></div>')
    elif status == "skipped":
        note = html_lib.escape(panel.get("note", "not available"))
        body = (f'<div class="state-card"><strong>Section not available</strong>'
                f'<p>{note}</p></div>')
    else:
        body = (f'<div class="exports">'
                f'<button onclick="reportExport(\'{pid}\',\'svg\')">Save SVG</button>'
                f'<button onclick="reportExport(\'{pid}\',\'png\')">Save PNG</button></div>'
                f'<div class="panel" id="panel-{pid}"></div>')
        body = f'<div class="figure-stage">{body}</div>'
    caption = html_lib.escape(panel.get("caption", ""))
    hidden = "" if active else " hidden"
    return (f'<section id="sec-{pid}" class="report-section" role="tabpanel" '
            f'aria-labelledby="tab-{pid}"{hidden}>\n'
            f'<div class="section-head"><h2>{title}</h2>{state}</div>\n{body}\n'
            f'<p class="caption">{caption}</p>\n</section>')


_RUNTIME = """
window.reportExport = function (panelId, kind) {
  const svg = document.querySelector("#panel-" + panelId + " svg");
  if (!svg) return;
  const base = document.querySelector("#panel-" + panelId + " canvas");
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  if (base) {
    const image = document.createElementNS("http://www.w3.org/2000/svg", "image");
    image.setAttribute("x", "0"); image.setAttribute("y", "0");
    image.setAttribute("width", svg.viewBox.baseVal.width || 1054);
    image.setAttribute("height", svg.viewBox.baseVal.height || 640);
    image.setAttribute("href", base.toDataURL("image/png"));
    clone.insertBefore(image, clone.firstChild);
  }
  const meta = (window.DATAByPanel[panelId] || {}).meta || {};
  const stamp = document.createComment(
    " chart-forge report | " + (meta.title || panelId) + " | " + (meta.source || "no-source") + " ");
  clone.insertBefore(stamp, clone.firstChild);
  const xml = '<?xml version="1.0" encoding="UTF-8"?>\\n' + new XMLSerializer().serializeToString(clone);
  const safe = s => (s || panelId).replace(/[^\\w.-]+/g, "_").slice(0, 80);
  const download = (name, href) => {
    const a = document.createElement("a");
    a.href = href; a.download = name; document.body.appendChild(a); a.click(); a.remove();
  };
  const url = URL.createObjectURL(new Blob([xml], { type: "image/svg+xml" }));
  if (kind === "svg") {
    download(safe(meta.title) + ".svg", url);
    setTimeout(() => URL.revokeObjectURL(url), 5000);
    return;
  }
  const img = new Image();
  img.onload = () => {
    const vb = svg.viewBox.baseVal;
    const w = vb.width || svg.clientWidth || 1054;
    const h = vb.height || svg.clientHeight || 600;
    const scale = 3.4;
    const canvas = document.createElement("canvas");
    canvas.width = Math.round(w * scale); canvas.height = Math.round(h * scale);
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, canvas.width, canvas.height);
    if (base) ctx.drawImage(base, 0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    URL.revokeObjectURL(url);
    download(safe(meta.title) + ".png", canvas.toDataURL("image/png"));
  };
  img.src = url;
};
window.reportSelect = function (panelId, detail) {
  const target = document.querySelector("#selection");
  if (!target) return;
  if (detail.kind === "family") {
    const ids = (detail.featureIds || []).slice(0, 8).join(", ");
    target.textContent = "Selected family " + (detail.familyId || detail.label || "") +
      (ids ? " · features: " + ids : "");
    return;
  }
  if (detail.kind === "sample") {
    target.textContent = "Selected sample " + (detail.sample || "") +
      (detail.group ? " · group: " + detail.group : "") +
      (Number.isFinite(Number(detail.x)) ? " · PC1 " + Number(detail.x).toFixed(3) : "") +
      (Number.isFinite(Number(detail.y)) ? " · PC2 " + Number(detail.y).toFixed(3) : "");
    return;
  }
  if (detail.kind === "feature") {
    target.textContent = "Selected feature " + (detail.featureId || "") +
      (Number.isFinite(Number(detail.log2fc)) ? " · log2FC " + Number(detail.log2fc).toFixed(3) : "") +
      (Number.isFinite(Number(detail.fdr)) ? " · BH FDR " + Number(detail.fdr).toExponential(2) : "");
    return;
  }
  const ids = (detail.assignmentIds || []).slice(0, 8).join(", ");
  const suffix = ids ? " · assignments: " + ids : "";
  target.textContent = "Selected " + detail.kind + " " + Number(detail.mass).toFixed(6) +
    " Da · " + detail.count + " PrSM observations" + suffix;
};
window.reportTab = function (panelId, updateHash) {
  const buttons = Array.from(document.querySelectorAll("[data-report-tab]"));
  const sections = Array.from(document.querySelectorAll(".report-section"));
  const target = "sec-" + panelId;
  let found = false;
  buttons.forEach(button => {
    const selected = button.dataset.reportTab === panelId;
    button.classList.toggle("is-active", selected);
    button.setAttribute("aria-selected", selected ? "true" : "false");
    if (selected) found = true;
  });
  sections.forEach(section => {
    const selected = section.id === target;
    section.hidden = !selected;
    section.setAttribute("aria-hidden", selected ? "false" : "true");
  });
  if (found && updateHash) history.replaceState(null, "", "#" + target);
};
(() => {
  const buttons = Array.from(document.querySelectorAll("[data-report-tab]"));
  if (!buttons.length) return;
  const requested = location.hash.replace(/^#sec-/, "");
  const initial = buttons.some(b => b.dataset.reportTab === requested)
    ? requested : buttons[0].dataset.reportTab;
  window.reportTab(initial, false);
  buttons.forEach((button, index) => {
    button.addEventListener("click", () => window.reportTab(button.dataset.reportTab, true));
    button.addEventListener("keydown", event => {
      if (event.key !== "ArrowRight" && event.key !== "ArrowLeft") return;
      event.preventDefault();
      const next = (index + (event.key === "ArrowRight" ? 1 : -1) + buttons.length) % buttons.length;
      buttons[next].focus();
      window.reportTab(buttons[next].dataset.reportTab, true);
    });
  });
})();
window.__chartReady = Boolean(document.querySelector(".panel svg, .panel canvas"));
"""

_REPORT_CSS = """
html { scroll-behavior: smooth; }
body { font-family: var(--font, system-ui); color: var(--ink, #1a1a1a);
       background: var(--c-gap, #f2f2f2); margin: 0; padding: 24px; }
main.report-shell { max-width: 1320px; min-height: calc(100vh - 48px); margin: 0 auto;
                    background: #fff; padding: 0 34px; box-shadow: 0 10px 36px rgba(26,26,26,.06); }
.masthead { padding: 30px 0 0; }
.eyebrow { margin: 0 0 8px; color: var(--c-link, #0072b2); font: 700 10px var(--font-mono, monospace);
           letter-spacing: .12em; text-transform: uppercase; }
h1 { font-size: clamp(1.65rem, 3vw, 2.35rem); line-height: 1.05; font-weight: 750;
     letter-spacing: -.045em; margin: 0 0 18px; }
dl.meta { display: grid; grid-template-columns: max-content minmax(0, 1fr) max-content minmax(0, 1fr);
          gap: .35rem 1rem; font-size: .76rem; color: var(--g-mid, #555); margin: 0;
          padding: 12px 0 16px; border-top: 1px solid var(--c-grid, #e0e0e0); }
dl.meta dt { font-weight: 700; color: var(--g-dark, #444); }
dl.meta dd { margin: 0; overflow-wrap: anywhere; }
.report-status { margin: 0; padding: 9px 0; font: 600 10px var(--font-mono, monospace);
                 color: var(--g-soft, #777); border-top: 1px solid var(--c-grid, #e0e0e0); }
nav.toc { position: sticky; top: 0; z-index: 5; display: flex; gap: 5px; overflow-x: auto;
          margin: 0 -34px; padding: 9px 34px; background: rgba(255,255,255,.97);
          border-top: 1px solid var(--c-grid, #e0e0e0); border-bottom: 1px solid var(--c-grid, #e0e0e0);
          scrollbar-width: thin; }
nav.toc button { flex: 0 0 auto; display: inline-flex; align-items: center; gap: 7px;
                 min-height: 34px; padding: 6px 12px; border: 1px solid transparent;
                 border-radius: 4px; background: transparent; color: var(--g-mid, #555);
                 font: 650 11px var(--font, system-ui); cursor: pointer; white-space: nowrap; }
nav.toc button:hover { color: var(--c-link, #0072b2); background: var(--c-gap, #f2f2f2); }
nav.toc button:focus-visible { outline: 2px solid var(--c-link, #0072b2); outline-offset: 2px; }
nav.toc button.is-active { color: var(--ink, #1a1a1a); border-color: var(--c-link, #0072b2);
                           background: #fff; }
nav.toc button.tab-primary { color: var(--c-link, #0072b2); }
nav.toc button.tab-primary.is-active { color: #fff; border-color: var(--c-link, #0072b2);
                                      background: var(--c-link, #0072b2); }
.tab-state { font-size: 9px; font-weight: 500; color: var(--g-faint, #999); }
nav.toc button.is-active .tab-state { color: var(--c-link, #0072b2); }
nav.toc button.tab-primary.is-active .tab-state { color: rgba(255,255,255,.82); }
#selection { min-height: 1.3em; margin: 10px 0 0; font-size: .76rem; color: var(--g-dark, #444);
             font-family: var(--font-mono, monospace); }
.report-section { margin: 0; padding: 30px 0 26px; scroll-margin-top: 58px; }
.section-head { display: flex; align-items: baseline; justify-content: space-between; gap: 16px;
                margin-bottom: 10px; }
.section-head h2 { font-size: 1.18rem; letter-spacing: -.02em; margin: 0; }
.state { flex: 0 0 auto; padding: 3px 8px; border-radius: 999px; font-size: 9px; font-weight: 650;
         letter-spacing: .02em; text-transform: uppercase; }
.state-ok { color: #176b4d; background: #e7f4ee; }
.state-skipped { color: #765c16; background: #fbf2d7; }
.state-error { color: #9c2e24; background: #fbe9e7; }
.exports { float: right; }
.exports button { font: 600 10.5px var(--font, system-ui); padding: 4px 9px;
                 border: 1px solid var(--c-grid, #e0e0e0); background: #fff;
                 color: var(--g-mid, #555); border-radius: 4px; cursor: pointer; }
.exports button:hover { color: var(--c-link, #0072b2); border-color: var(--c-link, #0072b2); }
.exports button:focus-visible { outline: 2px solid var(--c-link, #0072b2); outline-offset: 2px; }
.figure-stage { clear: both; margin-top: 12px; padding: 18px 14px 12px; background: #fff;
                border-top: 2px solid var(--ink, #1a1a1a); border-bottom: 1px solid var(--c-grid, #e0e0e0);
                overflow-x: auto; }
.panel { clear: both; min-width: 760px; text-align: center; }
.panel svg { max-width: 100%; height: auto; }
.caption { clear: both; font-size: .78rem; line-height: 1.5; color: var(--g-soft, #777); max-width: 70em;
           margin: 10px 0 0; }
.state-card { max-width: 680px; padding: 16px 18px; border: 1px solid #eadba9; background: #fffaf0;
              color: var(--g-dark, #444); }
.state-card strong { font-size: .86rem; }
.state-card p { margin: 6px 0 0; font-size: .78rem; color: var(--g-soft, #777); }
.state-card-error { border-color: #efc0ba; background: #fff7f6; }
footer { border-top: 1px solid var(--c-grid, #e0e0e0); margin-top: 18px;
         padding: 12px 0 30px; font-size: .72rem; color: var(--g-faint, #999); }
@media (max-width: 680px) {
  body { padding: 0; }
  main.report-shell { min-height: 100vh; padding: 0 16px; }
  nav.toc { margin-left: -16px; margin-right: -16px; padding-left: 16px; padding-right: 16px; }
  dl.meta { grid-template-columns: max-content minmax(0, 1fr); }
  .figure-stage { margin-left: -4px; margin-right: -4px; padding-left: 8px; padding-right: 8px; }
  .section-head { align-items: flex-start; }
  .section-head h2 { max-width: 74%; }
}
"""


def compile_report(panels, out_path, *, title: str, meta_rows=()) -> Path:
    """Assemble one self-contained report HTML.

    ``panels``: iterable of dicts. A chart panel needs
    ``{id, title, caption, script, data, meta?}`` where ``script`` is a Path
    whose code mounts via ``d3.select("#chart")`` and reads ``window.DATA``;
    ``data``/``meta`` become that panel's ``window.DATA``. A panel with
    ``status`` ("skipped"/"error") renders an in-page note instead — the human
    sees why a section is absent, not a silent hole.
    """
    panels = list(panels)
    data_by_panel = {}
    mounted_by_panel = {}
    for panel in panels:
        pid = panel["id"]
        if not re.fullmatch(r"[a-z0-9][a-z0-9-]*", pid):
            raise ValueError(f"panel id {pid!r} must be [a-z0-9-]")
        if panel.get("status"):
            continue
        script = _read(panel["script"])
        if MOUNT_TOKEN not in script:
            raise ValueError(f"panel {pid!r} script {panel['script']} never mounts "
                             f"{MOUNT_TOKEN}; the report cannot place it")
        mounted_by_panel[pid] = _js_safe(script.replace(MOUNT_TOKEN, f'd3.select("#panel-{pid}")'))
        data_by_panel[pid] = {"data": panel.get("data"),
                              "meta": {"title": panel["title"], **(panel.get("meta") or {})}}

    active_id = next((p["id"] for p in panels if p.get("status", "ok") == "ok"),
                     panels[0]["id"])
    body_panels = [
        _section(panel, active=panel["id"] == active_id) for panel in panels
    ]
    state_labels = {"ok": "ready", "skipped": "off", "error": "failed"}
    primary_ids = {"pca", "volcano"}
    toc = "".join(
        f'<button type="button" role="tab" class="{"tab-primary" if p["id"] in primary_ids else ""}" '
        f'id="tab-{html_lib.escape(p["id"], quote=True)}" '
        f'data-report-tab="{html_lib.escape(p["id"], quote=True)}" '
        f'aria-controls="sec-{html_lib.escape(p["id"], quote=True)}" '
        f'aria-selected="{"true" if p["id"] == active_id else "false"}">'
        f'{html_lib.escape(p["title"])} '
        f'<span class="tab-state">{html_lib.escape(state_labels.get(p.get("status", "ok"), "state"))}</span></button>'
        for p in panels
    )
    status_counts = {
        state: sum(1 for p in panels if p.get("status", "ok") == state)
        for state in ("ok", "skipped", "error")
    }
    status_text = (f'{len(panels)} views · {status_counts["ok"]} ready · '
                   f'{status_counts["skipped"]} unavailable · {status_counts["error"]} failed')
    meta_html = "".join(f"<dt>{html_lib.escape(k)}</dt><dd>{html_lib.escape(str(v))}</dd>"
                        for k, v in meta_rows)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    parts = [f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html_lib.escape(title)}</title>
<style>
{_read(FORGE / "templates" / "tokens.css")}
{_REPORT_CSS}
</style>
</head>
<body>
<main class="report-shell">
<header class="masthead">
<p class="eyebrow">TopPICR / run quality control</p>
<h1>{html_lib.escape(title)}</h1>
<dl class="meta">{meta_html}</dl>
<nav class="toc" role="tablist" aria-label="Report sections">{toc}</nav>
<p class="report-status">{html_lib.escape(status_text)}</p>
<p id="selection" aria-live="polite"></p>
</header>
"""]
    parts.extend(body_panels)
    parts.append(f'<footer>{html_lib.escape(title)} · generated {generated} · '
                 f'offline self-contained (vendored d3, no CDN)</footer>\n</main>\n')
    parts.append(f"<script>\n{_d3_library()}\n</script>\n")
    parts.append(f"<script>\n{_js_safe(_read(FORGE / 'templates' / 'micro.js'))}\n</script>\n")
    for provider in PROVIDER_SCRIPTS:
        parts.append(f"<script>\n{_provider(FORGE / provider)}\n</script>\n")
    parts.append(f"<script>\nwindow.DATAByPanel = {json_for_html(data_by_panel)};\n</script>\n")
    for panel in panels:
        if panel.get("status"):
            continue
        pid = panel["id"]
        parts.append(f"<script>\nwindow.DATA = window.DATAByPanel[{json_for_html(pid)}];\n</script>\n")
        parts.append(f"<script>\n{mounted_by_panel[pid]}\n</script>\n")
    parts.append(f"<script>\n{_RUNTIME}\n</script>\n</body>\n</html>\n")

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(parts), encoding="utf-8")
    return out


def envelope(report_path, panels) -> dict:
    """The agent-facing stdout contract shared by every report command:
    one JSON object with the report path and per-section status."""
    return {"report": str(report_path),
            "figures": sum(1 for p in panels if not p.get("status")),
            "sections": [{"id": p["id"], "status": p.get("status", "ok"),
                          **({"note": p["note"]} if p.get("note") else {})}
                         for p in panels]}
