#!/usr/bin/env python3
"""Absorb snow-d3 example HTMLs into self-contained offline versions.

The upstream examples import d3 from jsDelivr ESM; the vendored UMD build sets
window.d3 — so stripping the import line and inlining the vendored bundle makes
each example offline-loadable and pipeline-consistent (no CDN in scientific assets).

Usage: uv run scripts/absorb.py            # all examples -> assets/d3/examples-offline/
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

FORGE = Path(__file__).resolve().parent.parent
SRC = FORGE / "assets" / "d3" / "examples"
OUT = FORGE / "assets" / "d3" / "examples-offline"
VENDOR = FORGE / "assets" / "d3" / "vendor"
VENDORED = (VENDOR / "d3.min.js").read_text(encoding="utf-8")

IMPORT_RE = re.compile(
    r'import \* as d3 from "https://cdn\.jsdelivr\.net/npm/d3@[^"]*";?')


def _fix_sankey(html: str) -> str:
    """17-sankey: d3-sankey UMD extends global d3; shim the removed ESM namespace."""
    umd = (VENDOR / "d3-sankey.min.js").read_text(encoding="utf-8")
    html = html.replace(
        "<script>\n" + VENDORED + "\n</script>",
        "<script>\n" + VENDORED + "\n" + umd + "\n</script>\n<script>const sankey = d3;</script>")
    html = re.sub(r'import \* as sankey from "https://cdn\.jsdelivr\.net/npm/d3-sankey@[^"]*";?\n?', "", html)
    return html


def _fix_choropleth(html: str) -> str:
    """10-choropleth: embed world atlas JSON + topojson-client UMD, drop fetch/import."""
    topo = (VENDOR / "topojson-client.min.js").read_text(encoding="utf-8")
    world = (VENDOR / "countries-110m.json").read_text(encoding="utf-8")
    html = html.replace(
        "<script>\n" + VENDORED + "\n</script>",
        "<script>\n" + VENDORED + "\n" + topo + "\n</script>")
    html = html.replace(
        'const GEOJSON_URL = "https://cdn.jsdelivr.net/npm/world-atlas@2/countries-110m.json";',
        'const WORLD_110M = ' + world + ';')
    html = html.replace(
        """    const [topoMod, worldData] = await Promise.all([
      import("https://cdn.jsdelivr.net/npm/topojson-client@3/+esm"),
      fetch(GEOJSON_URL).then(r => r.json())
    ]);
    topojson = topoMod;
    topoData = worldData;""",
        """    topojson = window.topojson;
    topoData = WORLD_110M;""")
    return html


SPECIAL = {"17-sankey.html": _fix_sankey, "10-choropleth-map.html": _fix_choropleth}


def absorb(path: Path) -> str:
    html = path.read_text(encoding="utf-8")
    matches = IMPORT_RE.findall(html)
    if not matches:
        return f"SKIP {path.name} (no CDN import found)"
    html = IMPORT_RE.sub("", html)
    banner = ("<!-- absorbed by chart-forge/scripts/absorb.py — vendored d3 inlined, "
              "upstream: jiannanya/snow-d3 (MIT) -->\n")
    html = html.replace("<head>", "<head>\n" + banner, 1)
    # inject the vendored bundle right before the first module script
    html = html.replace("<script type=\"module\">",
                        "<script>\n" + VENDORED + "\n</script>\n<script type=\"module\">", 1)
    if path.name in SPECIAL:
        html = SPECIAL[path.name](html)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / path.name).write_text(html, encoding="utf-8")
    return f"OK {path.name} ({(OUT / path.name).stat().st_size / 1024:.0f} KB)"

if __name__ == "__main__":
    results = [absorb(p) for p in sorted(SRC.glob("*.html"))]
    print("\n".join(results))
    leftover = [p.name for p in OUT.glob("*.html")
                if "cdn.jsdelivr" in p.read_text(encoding="utf-8")]
    print(f"\n{sum(r.startswith('OK') for r in results)} absorbed, "
          f"{len(leftover)} with leftover CDN refs")
    sys.exit(1 if leftover else 0)
