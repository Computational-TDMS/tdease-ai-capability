// Report panel: scatter with optional series (colour+shape), line overlays
// (LOESS/fit curves) and dashed threshold reference lines.
// Contract: window.DATA.data = {
//   series: [{label, points: [{x, y, hover?}], shape?}],   // shape cycled if omitted
//   lines?: [{xs, ys, dash?}], vlines?: [{v, label?}], hlines?: [{v, label?}],
//   xlabel?, ylabel?, legend?: bool, height?: px, pointR?: px }
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, gridH } = forge;
  const { series, lines = [], vlines = [], hlines = [], xlabel, ylabel,
          legend = true, height, pointR = 3.4 } = window.DATA.data;
  const F = forge.fig({ W: 1054, H: height || 400 });
  const { W, H, k } = F;
  const M = { l: 64 * k, r: 24 * k, t: 26 * k, b: 52 * k };

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip();

  const SHAPES = [d3.symbolCircle, d3.symbolTriangle, d3.symbolDiamond,
                  d3.symbolSquare, d3.symbolCross];
  const all = series.flatMap(s => s.points);
  const xs = [...all.map(p => p.x), ...lines.flatMap(l => l.xs),
              ...vlines.map(l => l.v), ...hlines.map(l => l.v)];
  const ys = [...all.map(p => p.y), ...lines.flatMap(l => l.ys),
              ...vlines.map(l => l.v), ...hlines.map(l => l.v)];
  const [xlo, xhi] = d3.extent(xs), [ylo, yhi] = d3.extent(ys);
  const pad = (lo, hi2) => (hi2 - lo) * 0.05 || Math.abs(lo) * 0.05 || 1;
  const px = pad(xlo, xhi), py = pad(ylo, yhi);
  const x = d3.scaleLinear().domain([xlo - px, xhi + px]).range([M.l, W - M.r]);
  const y = d3.scaleLinear().domain([ylo - py, yhi + py]).range([H - M.b, M.t]);
  const fmt = d3.format(",.4~g");

  gridH(svg.append("g"), y, { w: W - M.l - M.r, l: M.l, n: 6 });

  for (const l of vlines) {
    svg.append("line").attr("x1", x(l.v)).attr("x2", x(l.v))
      .attr("y1", M.t).attr("y2", y.range()[0])
      .attr("stroke", C.threshold).attr("stroke-width", 1).attr("stroke-dasharray", "5 4");
    if (l.label) svg.append("text").attr("x", x(l.v)).attr("y", M.t - 6)
      .attr("text-anchor", "middle").attr("font-size", 10 * k).attr("fill", C.gMid).text(l.label);
  }
  for (const l of hlines) {
    svg.append("line").attr("x1", M.l).attr("x2", x.range()[1])
      .attr("y1", y(l.v)).attr("y2", y(l.v))
      .attr("stroke", C.threshold).attr("stroke-width", 1).attr("stroke-dasharray", "5 4");
    if (l.label) svg.append("text").attr("x", x.range()[1] - 4).attr("y", y(l.v) - 5)
      .attr("text-anchor", "end").attr("font-size", 10 * k).attr("fill", C.gMid).text(l.label);
  }

  for (const line of lines) {
    svg.append("path").datum(d3.zip(line.xs, line.ys))
      .attr("fill", "none").attr("stroke", C.gDark).attr("stroke-width", 1.6)
      .attr("stroke-dasharray", line.dash ? "6 4" : null)
      .attr("d", d3.line().x(d => x(d[0])).y(d => y(d[1])));
  }

  series.forEach((s, i) => {
    const color = s.color || C.cat[i % C.cat.length];
    const shape = SHAPES[(s.shapeIndex ?? i) % SHAPES.length];
    svg.append("g").selectAll("path").data(s.points).join("path")
      .attr("transform", p => `translate(${x(p.x)},${y(p.y)})`)
      .attr("d", d3.symbol().type(shape).size(Math.PI * pointR * pointR))
      .attr("fill", color).attr("fill-opacity", 0.85)
      .attr("stroke", "#fff").attr("stroke-width", 0.5)
      .style("cursor", "crosshair")
      .on("mousemove", (e, p) => showTip(tt, e, p.hover ||
        `<b>${s.label}</b><br>x ${fmt(p.x)} · y ${fmt(p.y)}`))
      .on("mouseleave", () => hideTip(tt));
  });

  axBottom(svg.append("g"), x, { y: y.range()[0], label: xlabel });
  axLeft(svg.append("g"), y, { l: M.l, label: ylabel });

  if (legend && series.length > 1) {
    const lg = svg.append("g").attr("transform",
      `translate(${W - M.r - 128 * k},${M.t})`);
    series.forEach((s, i) => {
      const color = s.color || C.cat[i % C.cat.length];
      const shape = SHAPES[(s.shapeIndex ?? i) % SHAPES.length];
      const g = lg.append("g").attr("transform", `translate(0,${i * 15 * k})`);
      g.append("path").attr("transform", "translate(6,0)")
        .attr("d", d3.symbol().type(shape).size(30))
        .attr("fill", color).attr("stroke", "#fff").attr("stroke-width", 0.5);
      g.append("text").attr("x", 18).attr("y", 4).attr("font-size", 10.5 * k)
        .attr("fill", C.gDark).text(s.label);
    });
  }
})();
