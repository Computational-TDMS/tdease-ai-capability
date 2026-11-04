;(function () {
  "use strict";
  const data = window.DATA.data || {};
  const points = data.points || [];
  const C = forge.C;
  const { W, H } = forge.fig({ W: 900, H: 620 });
  const M = { t: 52, r: 190, b: 66, l: 78 };
  const plotW = W - M.l - M.r, plotH = H - M.t - M.b;
  const svg = d3.select("#chart").append("svg")
    .attr("viewBox", `0 0 ${W} ${H}`).attr("width", W).attr("height", H)
    .style("font-family", C.font);
  const extent = values => {
    const e = d3.extent(values.filter(Number.isFinite));
    return e[0] === e[1] ? [e[0] - 1, e[1] + 1] : e;
  };
  const x = d3.scaleLinear().domain(extent(points.map(d => +d.x))).nice()
    .range([M.l, M.l + plotW]);
  const y = d3.scaleLinear().domain(extent(points.map(d => +d.y))).nice()
    .range([M.t + plotH, M.t]);
  const groups = Array.from(new Set(points.map(d => d.group || "other")));
  const color = d3.scaleOrdinal(groups, groups.map((_, i) => C.cat[i % C.cat.length]));
  const esc = value => String(value ?? "").replace(/[&<>\"]/g, c => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;"
  }[c]));
  const tip = forge.tip();
  const sizeValues = points.map(d => +d.sizeValue).filter(Number.isFinite);
  const sizeExtent = d3.extent(sizeValues);
  const sizeScale = sizeValues.length && sizeExtent[0] !== sizeExtent[1]
    ? d3.scaleSqrt().domain(sizeExtent).range([42, 230]) : null;
  const radius = d => Math.sqrt((sizeScale ? sizeScale(+d.sizeValue) : 150) / Math.PI);

  forge.gridV(svg, x, { h: plotH, t: M.t, n: 6 });
  forge.gridH(svg, y, { w: plotW, l: M.l, n: 6 });
  [0].forEach(v => {
    if (x.domain()[0] <= v && v <= x.domain()[1]) svg.append("line")
      .attr("x1", x(v)).attr("x2", x(v)).attr("y1", M.t).attr("y2", M.t + plotH)
      .attr("stroke", C.threshold).attr("stroke-dasharray", "5 4");
    if (y.domain()[0] <= v && v <= y.domain()[1]) svg.append("line")
      .attr("x1", M.l).attr("x2", M.l + plotW).attr("y1", y(v)).attr("y2", y(v))
      .attr("stroke", C.threshold).attr("stroke-dasharray", "5 4");
  });

  groups.forEach(group => {
    const pts = points.filter(d => (d.group || "other") === group)
      .map(d => [x(+d.x), y(+d.y)]);
    const hull = pts.length >= 3 ? d3.polygonHull(pts) : null;
    if (hull) svg.append("path")
      .attr("d", d3.line().curve(d3.curveCatmullRomClosed)(hull))
      .attr("fill", color(group)).attr("fill-opacity", 0.07)
      .attr("stroke", color(group)).attr("stroke-opacity", 0.35)
      .attr("stroke-dasharray", "4 3");
  });

  const marks = svg.append("g").selectAll("circle").data(points).join("circle")
    .attr("cx", d => x(+d.x)).attr("cy", d => y(+d.y)).attr("r", radius)
    .attr("fill", d => color(d.group || "other")).attr("fill-opacity", 0.82)
    .attr("stroke", "#fff").attr("stroke-width", 1.4)
    .style("cursor", "pointer")
    .on("mouseenter", function (event, d) {
      d3.select(this).attr("stroke", C.ink).attr("stroke-width", 2.2);
      const detail = d.hover || `<b>${esc(d.sample)}</b><br>${esc(d.group || "other")}`;
      forge.showTip(tip, event, detail, 300);
    })
    .on("mousemove", (event) => forge.showTip(tip, event, tip.html(), 300))
    .on("mouseleave", function () {
      d3.select(this).attr("stroke", "#fff").attr("stroke-width", 1.4);
      forge.hideTip(tip);
    })
    .on("click", (_, d) => window.reportSelect("pca", {
      kind: "sample", sample: d.sample, group: d.group, x: d.x, y: d.y
    }));

  forge.axBottom(svg, x, { y: M.t + plotH, n: 6, label: data.xlabel || "PC1" });
  forge.axLeft(svg, y, { l: M.l, n: 6, label: data.ylabel || "PC2" });

  const legend = svg.append("g").attr("transform", `translate(${M.l + plotW + 22},${M.t})`);
  legend.append("text").attr("font-size", 11).attr("font-weight", 700)
    .attr("fill", C.ink).text("Group");
  groups.forEach((group, i) => {
    const row = legend.append("g").attr("transform", `translate(0,${22 + i * 24})`);
    row.append("circle").attr("cx", 5).attr("cy", 0).attr("r", 5)
      .attr("fill", color(group)).attr("fill-opacity", 0.82);
    row.append("text").attr("x", 17).attr("y", 4).attr("font-size", 10)
      .attr("fill", C.gMid).text(`${group} (n=${points.filter(d => (d.group || "other") === group).length})`);
  });
  if (sizeScale) {
    const sy = 32 + groups.length * 24;
    legend.append("text").attr("y", sy + 12).attr("font-size", 11).attr("font-weight", 700)
      .attr("fill", C.ink).text(data.sizeLabel || "size");
    const refs = [sizeExtent[0], d3.median(sizeValues), sizeExtent[1]];
    refs.forEach((value, i) => {
      const r = radius({ sizeValue: value });
      const row = legend.append("g").attr("transform", `translate(${r + 4},${sy + 38 + i * 30})`);
      row.append("circle").attr("cy", -r).attr("r", r).attr("fill", C.gSoft).attr("fill-opacity", 0.22)
        .attr("stroke", C.gMid).attr("stroke-width", 0.7);
      row.append("text").attr("x", 18).attr("y", 4 - r).attr("font-size", 9)
        .attr("fill", C.gMid).text(d3.format(",.0f")(value));
    });
  }
  forge.hideTip(tip);
  window.__chartReady = true;
})();
