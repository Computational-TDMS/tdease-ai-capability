/* Family-10 fragment-annotation composite — the FIRST components-architecture
 * example: TopPIC PRSM title + annotated spectrum + ppm error panel, three
 * panels on one forge.panelStack sharing the m/z x-scale, one sync() choke
 * point for the three-way cleavage ↔ peak ↔ error-linkage. The title is a
 * static panel (zoom-invariant); the spectrum and error panels are wiped and
 * re-drawn by the stack on zoom. Distilled from the frontend
 * AnnotatedSpectrumIPSA.vue (IPSA/TopMSV lineage); fictional example data. */
(function () {
  "use strict";
  const { peaks, residues, tolerance_ppm: TOL, scan, annotations = [] } = window.DATA.data;
  const C = forge.C;
  const M = forgeMs2;
  const L = residues.length;
  const frags = peaks;                     // match table: every peak IS a fragment
  const basePeak = d3.max(frags, f => f.intensity);
  const rel = f => f.intensity / basePeak * 100;

  const F = forge.fig({ H: 860 });   // journal builds take preset width; k scales margins
  const k = F.k;
  const W = F.W, ML = Math.round(96 * k), MR = Math.round(40 * k);
  const specH = Math.round(420 * k), errH = Math.round(140 * k);
  const fs = { tick: F.J ? F.J.fsTick : 11, label: F.J ? F.J.fsLabel : 12, min: F.J ? F.J.fsMin : 10.5 };

  const root = d3.select("#chart").style("font-family", C.font);
  const svg = root.append("svg");          // sized after the stack measures its panels

  /* ---------- shared scales ---------- */
  const [m0, m1] = d3.extent(frags, f => f.mz);
  const pad = (m1 - m0) * 0.02;
  const x = d3.scaleLinear().domain([m0 - pad, m1 + pad]).range([ML, W - MR]);
  const fullX = x.domain();
  const y = d3.scaleLinear().domain([0, 105]).range([0, 100]);   // placeholder, set per draw

  /* ---------- controls (TopMSV peak-links idiom) ---------- */
  const controls = root.insert("div", "svg").style("display", "flex")
    .style("align-items", "center").style("gap", "14px")
    .style("margin", "0 0 6px 2px").style("font-size", "12px");
  let mode = "all";                        // all | matched | unmatched
  let showLabels = true;
  const links = {};
  [["all", "All"], ["matched", "Matched"], ["unmatched", "Not matched"]].forEach(([k, t]) => {
    links[k] = controls.append("span").style("cursor", "pointer")
      .style("color", C.link).style("font-weight", 700)
      .on("click", () => { mode = k; restyleLinks(); redraw(); });
  });
  function restyleLinks() {
    const counts = { all: frags.length, matched: frags.length, unmatched: 0 };
    Object.entries(links).forEach(([k, s]) => s
      .style("color", mode === k ? C.link : C.gMid)
      .style("font-weight", mode === k ? 700 : 400)
      .text(`${t[k]} (${counts[k]})`));
  }
  const t = { all: "All", matched: "Matched", unmatched: "Not matched" };
  const annLabel = controls.append("label").style("display", "inline-flex")
    .style("align-items", "center").style("gap", "5px").style("cursor", "pointer");
  annLabel.append("input").attr("type", "checkbox").property("checked", true)
    .on("change", e => { showLabels = e.target.checked; redraw(); });
  annLabel.append("span").text("annotations");
  const rangeText = controls.append("span").style("color", C.gSoft).style("font-size", "11px");
  svg.on("dblclick", reset);   // house convention: double-click = reset
  controls.append("span").style("color", C.gFaint).style("font-size", "10.5px")
    .text(`${scan} \u00b7 tolerance ${TOL} ppm`);

  /* ---------- tooltip ---------- */
  const tip = forge.tip();
  const fragTip = f =>
    `<b style="color:${M.ionColor(f.ion)}">${f.label}</b> \u00b7 charge ${f.charge}` +
    `<br>m/z ${f.mz.toFixed(4)} (theo ${f.theo.toFixed(4)})` +
    `<br>${f.ppm >= 0 ? "+" : ""}${f.ppm.toFixed(2)} ppm \u00b7 ${rel(f).toFixed(1)}% BPI`;

  /* ---------- panel bodies (the stack wipes each <g> before every draw) ---------- */
  let title, peakSel = null, curSw = 1.6, specTop = 0, errRef = null;

  function visible() {
    const [a, b] = x.domain();
    return frags.filter(f => f.mz >= a && f.mz <= b);  // match table: all matched
  }

  function drawSpectrum(g) {
    const vis = visible();
    const yMax = (vis.length ? d3.max(vis, rel) : 100) * 1.06;
    y.domain([0, yMax]);

    const gGrid = g.append("g");
    [0.25, 0.5, 0.75, 1].forEach(v => gGrid.append("line")
      .attr("x1", ML).attr("x2", W - MR)
      .attr("y1", y(v * yMax / 1.06)).attr("y2", y(v * yMax / 1.06))
      .attr("stroke", C.grid).attr("stroke-width", C.swGrid));

    M.annBands(g.append("g"), annotations, { x, top: specTop, height: specH, fs: fs.min });

    const gPeaks = g.append("g");
    const sw = Math.max(1, Math.min(2, (W - ML - MR) / Math.max(vis.length, 400) * 2.5));
    curSw = sw;
    peakSel = M.stems(gPeaks, vis, { x, y: f => y(rel(f)), y0: y(0), sw });
    M.peakHits(gPeaks, vis, {
      x, top: specTop, height: specH - 24,
      onMove: (e, f) => sync(frags.filter(g => g.ion === f.ion && g.number === f.number), e),
      onOut: () => sync(null),
    });

    drawLabels(g);

    const gAxes = g.append("g");
    gAxes.append("g").attr("transform", `translate(${ML},0)`)
      .call(d3.axisLeft(y).ticks(4).tickFormat(v => Math.round(v / yMax * 100) + "%").tickSize(3))
      .call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    gAxes.append("line").attr("x1", ML).attr("x2", W - MR)
      .attr("y1", y(0)).attr("y2", y(0)).attr("stroke", C.gFaint);
    gAxes.append("text").attr("x", ML - 44).attr("y", (specTop + specTop + specH - 24) / 2)
      .attr("transform", `rotate(-90 ${ML - 44} ${(specTop + specTop + specH - 24) / 2})`)
      .attr("text-anchor", "middle").attr("font-size", 11).attr("fill", C.gMid)
      .text("relative intensity (% BPI)");
  }

  function drawLabels(g) {
    if (showLabels) {
      M.topLabels(g.append("g"), visible(), { x, yTip: f => y(rel(f)) });
    }
    M.annCallouts(g.append("g"), annotations, { x, top: specTop, baseline: y(0), fs: fs.min });
  }

  /* ---------- the stack: title (static) → spectrum → ppm error ---------- */
  const stack = forge.panelStack(svg, [
    { static: true, gap: 16, draw: (g, top) => {
      title = M.prsmTitle(g, residues, frags, {
        top, k,
        onCleavage: (e, frag) => sync([frag], e),
        onCleavageLeave: () => sync(null),
        onResidue: (e, pos) => {
          const covering = frags.filter(f =>
            f.ion === "c" ? pos <= f.number : pos >= L - f.number + 1);
          sync(covering, e, pos);
        },
        onResidueLeave: () => sync(null),
      });
      return title.height;
    } },
    { gap: 26, draw: (g, top) => {
      specTop = top;
      y.range([specTop + specH - 24, specTop + 8]);
      drawSpectrum(g);
      return specH;
    } },
    { draw: (g, top) => {
      errRef = M.errorPanel(g, frags, { top, height: errH, x, tol: TOL });
      errRef.dots
        .on("mouseenter", (e, f) => sync([f], e))
        .on("mouseleave", () => sync(null));
      return errH;
    } },
  ], { top: Math.round(12 * k), W, k });

  svg.attr("viewBox", `0 0 ${W} ${stack.H + 8}`).attr("width", W).attr("height", stack.H + 8);

  /* ---------- three-way sync choke point ----------
   * sel: fragment array (or null); resPos: the residue that fired, if any.
   * Dims every panel consistently: selected cleavages/peaks/dots stay lit,
   * covered residues go accent. */
  const coverOf = f => f.ion === "c"
    ? d3.range(1, f.number + 1) : d3.range(L - f.number + 1, L + 1);
  function sync(sel, e, resPos) {
    const keys = sel ? new Set(sel.map(M.fkey)) : null;
    const covPos = sel ? new Set(sel.flatMap(coverOf)) : null;

    title.marks
      .attr("stroke-width", k => keys && keys.has(k) ? 2.4 : 1.2)
      .attr("opacity", k => !keys ? 0.95 : (keys.has(k) ? 1 : 0.15));
    title.letters
      .attr("font-weight", d => (d.mod || (covPos && covPos.has(d.pos))) ? 700 : 500)
      .attr("fill", d => d.mod ? C.accent : (covPos && covPos.has(d.pos) ? C.link : C.ink));

    if (peakSel) peakSel
      .attr("opacity", f => (!keys ? 0.92 : (keys.has(M.fkey(f)) ? 1 : 0.18)))
      .attr("stroke-width", f => keys && keys.has(M.fkey(f)) ? curSw + 1.4 : curSw);
    if (errRef) errRef.dots
      .attr("r", f => keys && keys.has(M.fkey(f)) ? 4.2 : 2.6)
      .attr("opacity", f => !keys ? 0.85 : (keys.has(M.fkey(f)) ? 1 : 0.15));

    if (e) {
      let html;
      if (resPos) {
        const r = residues.find(q => q.pos === resPos);
        html = `<b>pos ${resPos} \u00b7 ${r.char}</b>` +
          (r.mod ? `<br><span style="color:${C.accent}">mod ${r.mod}</span>` : "") +
          `<br>${sel.length} fragment${sel.length === 1 ? "" : "s"} cover this residue`;
      } else {
        html = sel.map(fragTip).join("<br>\u2014<br>");
      }
      forge.showTip(tip, e, html, 300);
    } else {
      forge.hideTip(tip);
    }
  }

  /* ---------- brush zoom (micro controller; the stack re-draws panels) ---------- */
  forge.zoomBrush(svg, {
    extent: [[ML, specTop], [W - MR, specTop + specH - 24]],
    x,
    onZoom: dom => {
      x.domain(dom);
      window.__zoomState = {
        mz: dom.map(v => +v.toFixed(2)),
        zoomed: String(dom) !== String(fullX),
      };
      redraw();
    },
  });

  function redraw() {
    stack.redraw();
    rangeText.text(`m/z ${x.domain()[0].toFixed(1)} \u2013 ${x.domain()[1].toFixed(1)}`);
  }
  function reset() {
    x.domain(fullX);
    window.__zoomState = { mz: fullX.map(v => +v.toFixed(2)), zoomed: false };
    redraw();
  }

  restyleLinks();
  reset();
  window.__chartReady = true;
})();
