/* LC-MS feature map: interval bars + 2D brush zoom.
 * Contract: window.DATA.data = {features:[{mz, rt0, rt1, apex, logI, q,
 * dataset?, featureId?, mass?}], mode:"interval", n, sourceCount}.
 * A feature is drawn as a horizontal elution-window bar; its apex is a notch.
 */
(function () {
  "use strict";
  const payload = window.DATA.data || {};
  const feats = payload.features || [];
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
  const fullX = d3.scaleLinear()
    .domain([d3.min(feats, d => d.rt0), d3.max(feats, d => d.rt1)]).nice();
  const fullY = d3.scaleLinear().domain(d3.extent(feats, d => d.mz)).nice();
  const x = d3.scaleLinear().domain(fullX.domain()).range([M.l, W - M.r]);
  const y = d3.scaleLinear().domain(fullY.domain()).range([H - M.b, M.t]);
  const [iLo, iHi] = d3.extent(feats, d => d.logI);
  const ramp = d3.interpolateRgbBasis([...batlow].reverse());
  const color = d3.scaleSequential(ramp).domain([iLo, iHi]);
  const ordered = [...feats].sort((a, b) => a.logI - b.logI);
  const barH = 1.7, apexH = 3.4;

  function drawCloud() {
    ctx.clearRect(0, 0, W, H);
    ctx.globalAlpha = 0.9;
    for (const d of ordered) {
      const x0 = x(d.rt0), x1 = x(Math.max(d.rt1, d.rt0 + 0.05));
      const cy = y(d.mz);
      ctx.fillStyle = color(d.logI);
      ctx.fillRect(x0, cy - barH / 2, Math.max(0.7, x1 - x0), barH);
      ctx.fillRect(x(d.apex) - 0.75, cy - apexH / 2, 1.5, apexH);
    }
    ctx.globalAlpha = 1;
  }
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
      .call(d3.axisBottom(x).ticks(8)).call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    gAxes.append("g").attr("transform", `translate(${M.l},0)`)
      .call(d3.axisLeft(y).ticks(6)).call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    gAxes.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", H - 10)
      .attr("text-anchor", "middle").attr("font-size", 13).attr("fill", C.gMid).text("elution window (min)");
    gAxes.append("text").attr("transform", "rotate(-90)").attr("x", -(M.t + H - M.b) / 2)
      .attr("y", 16).attr("text-anchor", "middle").attr("font-size", 13).attr("fill", C.gMid).text("average m/z");
  }
  const cbW = 300, cbX = (W - cbW) / 2, cbY = 6;
  const defs = svg.append("defs");
  const gradient = defs.append("linearGradient").attr("id", "featuremap-intensity");
  d3.range(31).forEach(i => gradient.append("stop").attr("offset", i * 100 / 30 + "%")
    .attr("stop-color", color(iLo + (iHi - iLo) * i / 30)));
  svg.append("rect").attr("x", cbX).attr("y", cbY).attr("width", cbW).attr("height", 9)
    .attr("fill", "url(#featuremap-intensity)").attr("rx", 2);
  svg.append("text").attr("x", cbX - 8).attr("y", cbY + 8).attr("text-anchor", "end")
    .attr("font-size", 10).attr("fill", C.gMid).text(`10^${iLo.toFixed(1)}`);
  svg.append("text").attr("x", cbX + cbW + 8).attr("y", cbY + 8)
    .attr("font-size", 10).attr("fill", C.gMid).text(`10^${iHi.toFixed(1)}`);
  svg.append("text").attr("x", cbX + cbW / 2).attr("y", cbY + 22)
    .attr("text-anchor", "middle").attr("font-size", 10).attr("fill", C.gMid)
    .text("apex intensity · darker = stronger");
  const zoomLabel = svg.append("text").attr("x", W - M.r).attr("y", cbY + 8)
    .attr("text-anchor", "end").attr("font-size", 10).attr("fill", C.gMid);
  function applyZoom(xDomain, yDomain) {
    x.domain(xDomain); y.domain(yDomain); drawCloud(); renderGrid(); renderAxes();
    const zoomed = String(xDomain) !== String(fullX.domain()) || String(yDomain) !== String(fullY.domain());
    zoomLabel.text(zoomed ? `RT ${xDomain[0].toFixed(1)}–${xDomain[1].toFixed(1)} · m/z ${yDomain[0].toFixed(0)}–${yDomain[1].toFixed(0)}` : "");
    window.__zoomState = { rt: xDomain, mz: yDomain, zoomed };
  }
  const brush = d3.brush().extent([[M.l, M.t], [W - M.r, H - M.b]])
    .on("end", event => {
      if (!event.selection) return;
      const [[x0, y0], [x1, y1]] = event.selection;
      if (x1 - x0 < 8 || y1 - y0 < 8) { svg.select(".brush").call(brush.move, null); return; }
      applyZoom([x.invert(x0), x.invert(x1)], [y.invert(y1), y.invert(y0)]);
      svg.select(".brush").call(brush.move, null);
    });
  svg.append("g").attr("class", "brush").call(brush).selectAll(".selection")
    .attr("fill", C.link).attr("fill-opacity", 0.08).attr("stroke", C.link).attr("stroke-width", C.swAxis);
  svg.on("dblclick", () => applyZoom(fullX.domain(), fullY.domain()));

  const tooltip = forge.tip();
  svg.on("mousemove", (event) => {
    const [mx, my] = d3.pointer(event);
    if (mx < M.l || mx > W - M.r || my < M.t || my > H - M.b) { forge.hideTip(tooltip); return; }
    let best = null, bestD = Infinity;
    for (const d of feats) {
      const dx = Math.max(x(d.rt0) - mx, 0, mx - x(d.rt1));
      const dy = y(d.mz) - my;
      const dist = dx * dx + dy * dy;
      if (dist < bestD) { bestD = dist; best = d; }
    }
    if (!best || bestD > 900) { forge.hideTip(tooltip); return; }
    forge.showTip(tooltip, event,
      `<b>${best.dataset || "feature"}</b> · feature ${best.featureId ?? "—"}<br>` +
      `m/z ${Number(best.mz).toFixed(5)} · charge ${best.q ?? "—"}<br>` +
      `elution ${best.rt0}–${best.rt1} min · apex ${best.apex}<br>` +
      `mass ${best.mass == null ? "—" : Number(best.mass).toFixed(4)} Da · apex intensity ${(10 ** best.logI).toExponential(2)}`);
  }).on("mouseleave", () => forge.hideTip(tooltip));
  drawCloud(); renderGrid(); renderAxes(); applyZoom(fullX.domain(), fullY.domain());
  window.__chartReady = true;
})();
