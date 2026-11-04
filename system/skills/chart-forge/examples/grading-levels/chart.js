// Smith proteoform grading histogram — grouped bars per level 1-5, counting
// assignments as searched and after the pipeline's direct attribution.
// Statistics view: nothing on this chart is a chemical claim. Contract:
// { data: { levels: [{ level, searched, attributed }],
//           series: [{ key, label }] }, meta }
(function () {
  "use strict";
  const { C, tip, showTip, hideTip, axBottom, axLeft, gridH, legend } = forge;
  const levels = window.DATA.data.levels;
  const series = window.DATA.data.series || [
    { key: "searched", label: "as searched" }];
  const F = forge.fig({ W: 1054, H: 400 });
  const { W, H, k } = F;
  const M = { l: 64 * k, r: 42 * k, t: 52 * k, b: 56 * k };

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  const colorOf = { searched: C.gGhost, attributed: C.link };
  const yMax = d3.max(levels, d => d3.max(series, s => d[s.key]));
  const yMin = Math.min(0, d3.min(levels, d => d3.min(series, s => d[s.key])));
  const x0 = d3.scaleBand().domain(levels.map(d => String(d.level)))
    .range([M.l, W - M.r]).paddingInner(0.4);
  const x1 = d3.scaleBand().domain(series.map(s => s.key))
    .range([0, x0.bandwidth()]).padding(0.12);
  const y = d3.scaleLinear().domain([yMin * 1.15, yMax * 1.15]).range([H - M.b, M.t]);

  gridH(svg.append("g"), y, { w: W - M.l - M.r, l: M.l, n: 5 });
  const g = svg.append("g");
  levels.forEach(d => {
    series.forEach(s => {
      g.append("rect")
        .attr("x", x0(String(d.level)) + x1(s.key)).attr("y", y(d[s.key]))
        .attr("width", x1.bandwidth())
        .attr("height", Math.abs(y(0) - y(d[s.key])))
        .attr("fill", colorOf[s.key] || C.cat[series.indexOf(s) % C.cat.length])
        .attr("fill-opacity", 0.85)
        .style("cursor", "crosshair")
        .on("mousemove", e => showTip(tt, e,
          `<b>Level ${d.level}</b> · ${s.label}<br>${d[s.key]} assignments`, 220))
        .on("mouseleave", () => hideTip(tt));
    });
  });

  legend(svg, series.map(s => ({ label: s.label, color: colorOf[s.key] || C.cat[0] })),
         { x: M.l, y: 24 * k, pitch: 190 * k });
  axLeft(svg, y, { l: M.l, n: 5, label: "assignments" });
  axBottom(svg, x0, { y: H - M.b + 6 * k, label:
    "Smith proteoform level (1 = fully defined … 5 = observed mass only)" });
  window.__chartReady = true;
})();
