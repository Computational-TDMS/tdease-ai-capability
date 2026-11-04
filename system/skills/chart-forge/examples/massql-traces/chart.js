// massql plot.py, reproduced — MassQL SELECT traces on a shared RT axis.
// Two stacked magnitude domains (TIC ~1e8 vs XIC ~1e4) instead of the
// matplotlib single raw axis, so neither trace collapses the other.
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, fmt, zoomBrush } = forge;
  const traces = window.DATA.data.traces;
  const F = forge.fig({ H: 560 });
  const { W, H, k } = F;
  const M = { l: 64 * k, r: 16 * k, t: 34 * k };

  const H_TIC = 165 * k, GAP = 40 * k;
  const yTopTic = M.t, yBotTic = M.t + H_TIC;
  const yTopXic = yBotTic + GAP, yBotXic = H - 34 * k;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  const tic = traces.find(t => t.kind === "tic");
  const xics = traces.filter(t => t.kind === "xic");
  const XIC_COLOR = [C.link, C.gMid];

  const full = d3.extent(traces.flatMap(t => t.rt));
  const x = d3.scaleLinear().domain(full).range([M.l, W - M.r]);
  const yT = d3.scaleLinear().domain([0, d3.max(tic.y) * 1.05]).nice().range([yBotTic, yTopTic]);
  const yX = d3.scaleLinear().domain([0, d3.max(xics.flatMap(t => t.y)) * 1.08]).nice().range([yBotXic, yTopXic]);

  /* ---------- TIC context panel */
  const gT = svg.append("g");
  const areaT = d3.area().x((d, i) => x(tic.rt[i])).y0(yBotTic).y1(d => yT(d));
  const pathArea = gT.append("path").attr("fill", C.gFaint).attr("fill-opacity", 0.35);
  const pathLine = gT.append("path").attr("fill", "none").attr("stroke", C.gSoft).attr("stroke-width", 1.2);
  svg.append("text").attr("x", M.l).attr("y", yTopTic - 8 * k)
    .attr("font-size", 11 * k).attr("font-weight", 700).attr("fill", C.gMid).text("TIC · MS1");
  axLeft(gT, yT, { l: M.l, n: 4, label: "total ion current" });

  /* ---------- XIC panel */
  const gX = svg.append("g");
  axLeft(gX, yX, { l: M.l, n: 5, label: "matched intensity" });
  const pathX = {};
  xics.forEach(t => { pathX[t.label] = gX.append("path").attr("fill", "none").attr("stroke-width", 1.6); });
  const on = {}; xics.forEach(t => on[t.label] = true);

  /* ---------- bottom axes: one per panel, rebound on zoom */
  const axT = axBottom(gT, x, { y: yBotTic, n: 10 });
  const axX = axBottom(gX, x, { y: yBotXic, n: 10, label: "retention time (min)" });
  function reaxis() {
    [axT, axX].forEach(sel => {
      const y = +sel.attr("transform").split(",")[1].replace(")", "");
      sel.call(d3.axisBottom(x).ticks(10));
      sel.select(".domain").attr("stroke-width", C.swAxis);
      sel.selectAll("text").attr("font-size", 11);
    });
  }

  function draw() {
    pathArea.attr("d", areaT(tic.y));
    pathLine.attr("d", d3.line()(tic.y.map((v, i) => [x(tic.rt[i]), yT(v)])));
    xics.forEach((t, ti) => {
      pathX[t.label]
        .attr("d", d3.line()(t.y.map((v, i) => [x(t.rt[i]), yX(v)])))
        .attr("stroke", XIC_COLOR[ti])
        .attr("display", on[t.label] ? null : "none");
    });
  }
  draw();

  /* ---------- legend with toggles */
  const leg = svg.append("g");
  let lx = M.l;
  xics.forEach((t, ti) => {
    const g = leg.append("g").attr("transform", `translate(${lx},${yTopXic - 16 * k})`).style("cursor", "pointer");
    g.append("line").attr("x2", 18 * k).attr("stroke", XIC_COLOR[ti]).attr("stroke-width", 2.6);
    g.append("text").attr("x", 24 * k).attr("y", 4 * k).attr("font-size", 11.5 * k)
      .attr("fill", C.gDark).attr("font-weight", 600).text(t.label);
    g.on("click", () => { on[t.label] = !on[t.label]; g.attr("opacity", on[t.label] ? 1 : 0.35); draw(); });
    lx += 24 * k + t.label.length * 6.8 * k + 30 * k;
  });

  /* ---------- hover: vertical rule + nearest values per trace */
  const bisect = d3.bisector((a, b) => a - b).left;
  const rule = svg.append("line").attr("y1", yTopTic).attr("y2", yBotXic)
    .attr("stroke", C.gGhost).attr("stroke-dasharray", "3 3").style("display", "none");
  svg.append("rect")
    .attr("x", M.l).attr("y", yTopTic).attr("width", W - M.r - M.l).attr("height", yBotXic - yTopTic)
    .attr("fill", "transparent")
    .on("mousemove", (e) => {
      const [mx] = d3.pointer(e);
      const rt = x.invert(mx);
      rule.attr("x1", mx).attr("x2", mx).style("display", null);
      const rows = [];
      traces.forEach(t => {
        if (t.kind === "xic" && !on[t.label]) return;
        const i = Math.min(Math.max(bisect(t.rt, rt), 0), t.rt.length - 1);
        rows.push(`<b>${t.label}</b><br>${fmt(t.y[i])} @ ${t.rt[i].toFixed(2)} min`);
      });
      showTip(tt, e, rows.join("<hr style='border:none;border-top:1px solid #ddd;margin:4px 0'>"), 250);
    })
    .on("mouseleave", () => { rule.style("display", "none"); hideTip(tt); });

  /* ---------- shared-RT brush zoom */
  const zoomTo = (dom) => { x.domain(dom); draw(); reaxis(); };
  zoomBrush(svg, {
    extent: [[M.l, yTopTic], [W - M.r, yBotXic]],
    x,
    onZoom: zoomTo,
  });
  svg.on("dblclick", () => zoomTo(full));   // house convention: double-click = reset
  window.__chartReady = true;
})();
