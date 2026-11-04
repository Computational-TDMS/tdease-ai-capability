/* LC-MS feature map v4 — elution-window bars + 2D brush zoom.
 * Each feature = horizontal bar [rt0, rt1] at its average m/z; colour = apex
 * intensity (batlow reversed). SVG d3.brush box-selects a region → canvas + axes
 * re-render to the selection; double-click or "Reset view" restores the full window.
 * Charge stays on hover only (representative value, not an intrinsic dimension).
 * Canvas hybrid renderer (>20k intervals). */
(function () {
  "use strict";
  const feats = window.DATA.data;
  const meta = window.DATA.meta;
  const C = forge.C;

  const holder = d3.select("#chart");
  const W = 1054, H = 640, M = { t: 34, r: 24, b: 52, l: 64 };
  holder.style("position", "relative").style("height", H + "px");
  const canvas = holder.append("canvas")
    .attr("width", W * 2).attr("height", H * 2)
    .style("width", W + "px").style("height", H + "px")
    .style("position", "absolute").style("left", 0).style("top", 0);
  const ctx = canvas.node().getContext("2d");
  ctx.scale(2, 2);
  const svg = holder.append("svg").attr("viewBox", `0 0 ${W} ${H}`)
    .attr("width", W).attr("height", H)
    .style("position", "absolute").style("left", 0).style("top", 0)
    .style("font-family", C.font);

  /* full (default) domains, nice — restored on reset */
  const fullX = d3.scaleLinear()
    .domain([d3.min(feats, d => d.rt0), d3.max(feats, d => d.rt1)]).nice();
  const fullY = d3.scaleLinear().domain(d3.extent(feats, d => d.mz)).nice();
  const x = d3.scaleLinear().domain(fullX.domain()).range([M.l, W - M.r]);
  const y = d3.scaleLinear().domain(fullY.domain()).range([H - M.b, M.t]);
  const [iLo, iHi] = d3.extent(feats, d => d.logI);
  const ramp = d3.interpolateRgbBasis([...batlow].reverse());
  const color = d3.scaleSequential(ramp).domain([iLo, iHi]);
  const colorLin = d3.scaleSequential(ramp).domain([10 ** iLo, 10 ** iHi]);
  const ordered = [...feats].sort((a, b) => a.logI - b.logI);

  /* ---- colour-scale toggle state (migrated from frontend FeatureMapViewer):
   * log (default) keeps weak features separable; linear shows the honest
   * dominance of the strong end. Same ramp, different domain.
   * Toggle UI is created after the colourbar frame (it needs cbY). ---- */
  let logScale = true;
  const cbStops = svg.append("defs").append("linearGradient").attr("id", "cb");
  function applyScale() {
    d3.range(31).forEach(i => cbStops.select(`stop:nth-child(${i + 1})`)
      .attr("stop-color", logScale
        ? color(iLo + (iHi - iLo) * i / 30)
        : colorLin(10 ** iLo + (10 ** iHi - 10 ** iLo) * i / 30)));
    drawCloud();
  }

  /* ---- canvas cloud (re-drawn on every zoom / scale toggle) ---- */
  const barH = 1.7, apexH = 3.4;
  function drawCloud() {
    ctx.clearRect(0, 0, W, H);
    ctx.globalAlpha = 0.9;
    for (const d of ordered) {
      const x0 = x(d.rt0), x1 = x(Math.max(d.rt1, d.rt0 + 0.05));
      const cy = y(d.mz);
      const col = logScale ? color(d.logI) : colorLin(10 ** d.logI);
      ctx.fillStyle = col;
      ctx.fillRect(x0, cy - barH / 2, x1 - x0, barH);
      ctx.fillRect(x(d.apex) - 0.75, cy - apexH / 2, 1.5, apexH);
    }
    ctx.globalAlpha = 1;
  }

  /* ---- SVG overlay: grid + axes re-render per zoom; colourbar static ---- */
  let gGrid = svg.append("g"), gAxes = svg.append("g");
  function renderGrid() {
    gGrid.remove(); gGrid = svg.insert("g", ":first-child");
    x.ticks(8).forEach(v => gGrid.append("line")
      .attr("x1", x(v)).attr("x2", x(v)).attr("y1", M.t).attr("y2", H - M.b)
      .attr("stroke", C.grid).attr("stroke-width", C.swGrid));
    y.ticks(6).forEach(v => gGrid.append("line")
      .attr("x1", M.l).attr("x2", W - M.r).attr("y1", y(v)).attr("y2", y(v))
      .attr("stroke", C.grid).attr("stroke-width", C.swGrid));
  }
  function renderAxes() {
    gAxes.remove(); gAxes = svg.append("g");
    gAxes.append("g").attr("transform", `translate(0,${H - M.b})`)
      .call(d3.axisBottom(x).ticks(8).tickSize(3))
      .call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    gAxes.append("g").attr("transform", `translate(${M.l},0)`)
      .call(d3.axisLeft(y).ticks(6).tickSize(3))
      .call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    gAxes.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", H - 10)
      .attr("text-anchor", "middle").attr("font-size", 13).text("elution window (min)");
    gAxes.append("text").attr("transform", "rotate(-90)").attr("x", -(M.t + H - M.b) / 2)
      .attr("y", 16).attr("text-anchor", "middle").attr("font-size", 13)
      .text("average m/z");
  }

  /* static colourbar frame (stops updated by applyScale) */
  const cbW = 300, cbX = (W - cbW) / 2, cbY = 6, cbH = 9;
  d3.range(31).forEach(i => cbStops.append("stop").attr("offset", i * (100 / 30) + "%"));
  const grad = cbStops;
  const cbLo = svg.append("text").attr("x", cbX - 8).attr("y", cbY + 8)
    .attr("text-anchor", "end").attr("font-size", 10).attr("fill", C.gMid);
  const cbHi = svg.append("text").attr("x", cbX + cbW + 8).attr("y", cbY + 8)
    .attr("font-size", 10).attr("fill", C.gMid);
  svg.append("rect").attr("x", cbX).attr("y", cbY).attr("width", cbW).attr("height", cbH)
    .attr("fill", "url(#cb)").attr("rx", 2);
  svg.append("text").attr("x", cbX + cbW / 2).attr("y", cbY + 22)
    .attr("text-anchor", "middle").attr("font-size", 10).attr("fill", C.gMid)
    .text("apex intensity (batlow reversed — darker = stronger)");

  /* toggle UI — after the colourbar frame (needs cbY); left of the bar */
  const toggle = svg.append("text").attr("x", M.l).attr("y", cbY + 8)
    .attr("font-size", 10).style("cursor", "pointer");
  const tLog = toggle.append("tspan").attr("fill", C.link).attr("font-weight", 700)
    .text("Log").on("click", () => { logScale = true; restyle(); });
  toggle.append("tspan").attr("fill", C.gFaint).text("  /  ");
  const tLin = toggle.append("tspan").attr("fill", C.gMid)
    .text("Linear").on("click", () => { logScale = false; restyle(); });
  function restyle() {
    tLog.attr("fill", logScale ? C.link : C.gMid)
      .attr("font-weight", logScale ? 700 : 400);
    tLin.attr("fill", !logScale ? C.link : C.gMid)
      .attr("font-weight", !logScale ? 700 : 400);
    applyScale();
  }
  function labelScale() {
    const fmt = v => logScale ? "10^" + v.toFixed(1) : v.toExponential(1);
    cbLo.text(fmt(logScale ? iLo : 10 ** iLo));
    cbHi.text(fmt(logScale ? iHi : 10 ** iHi));
  }

  /* zoom state indicator (hidden at full extent); double-click resets */
  const zoomLabel = svg.append("text").attr("x", W - M.r).attr("y", cbY + 8)
    .attr("text-anchor", "end").attr("font-size", 10).attr("fill", C.gMid);
  function applyZoom(xDomain, yDomain) {
    x.domain(xDomain); y.domain(yDomain);
    drawCloud(); renderGrid(); renderAxes();
    const zoomed = String(xDomain) !== String(fullX.domain()) ||
                   String(yDomain) !== String(fullY.domain());
    zoomLabel.text(zoomed
      ? `view: RT ${xDomain[0].toFixed(1)}–${xDomain[1].toFixed(1)} min · m/z ${yDomain[0].toFixed(0)}–${yDomain[1].toFixed(0)}`
      : "");
    /* deep-linkable view state (paradigm: interaction state must be shareable) */
    window.__zoomState = { rt: xDomain, mz: yDomain, zoomed };
  }

  /* ---- 2D brush: box-select → zoom ---- */
  const brush = d3.brush()
    .extent([[M.l, M.t], [W - M.r, H - M.b]])
    .on("end", (e) => {
      if (!e.selection) return;
      const [[x0, y0], [x1, y1]] = e.selection;
      if (x1 - x0 < 8 || y1 - y0 < 8) { svg.select(".brush").call(brush.move, null); return; }
      applyZoom([x.invert(x0), x.invert(x1)], [y.invert(y1), y.invert(y0)]);
      svg.select(".brush").call(brush.move, null);
    });
  svg.append("g").attr("class", "brush").call(brush)
    .selectAll(".selection")
    .attr("fill", C.link).attr("fill-opacity", 0.08)
    .attr("stroke", C.link).attr("stroke-width", C.swAxis);
  svg.on("dblclick", () => applyZoom(fullX.domain(), fullY.domain()));

  /* ---- hover: nearest feature by apex position; tooltip = full elution window ---- */
  const tooltip = forge.tip();
  svg.on("mousemove", (e) => {
    const [mx, my] = d3.pointer(e);
    if (mx < M.l || mx > W - M.r || my < M.t || my > H - M.b) { tooltip.style("opacity", 0); return; }
    let best = null, bestD = Infinity;
    for (const d of feats) {
      const dx = Math.max(x(d.rt0) - mx, 0, mx - x(d.rt1));
      const dy = y(d.mz) - my;
      const dist = dx * dx + dy * dy;
      if (dist < bestD) { bestD = dist; best = d; }
    }
    if (!best || bestD > 900) { tooltip.style("opacity", 0); return; }
    forge.showTip(tooltip, e,
      `m/z <b>${best.mz}</b> · charge ${best.q}<br>` +
      `elution ${best.rt0} – ${best.rt1} min (apex ${best.apex})<br>` +
      `apex intensity ≈ ${(10 ** best.logI).toExponential(1)}`);
  }).on("mouseout", () => tooltip.style("opacity", 0));

  drawCloud(); renderGrid(); renderAxes(); labelScale(); applyScale();
  window.__chartReady = true;
})();
