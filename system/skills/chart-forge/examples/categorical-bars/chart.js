/* Basic vertical bar chart — the minimal chart-forge starter. */
(function () {
  "use strict";
  const rows = window.DATA.data.rows;
  const C = forge.C;
  const F = forge.fig({ H: 320, double: true });
  const M = { l: 56, r: 20, t: 24, b: 42 };

  const svg = d3.select("#chart").append("svg")
    .attr("width", F.W).attr("height", F.H).attr("viewBox", `0 0 ${F.W} ${F.H}`);
  const tooltip = forge.tip();

  // Domain always includes 0, so bars have a real baseline on either side.
  const y = d3.scaleLinear()
    .domain([Math.min(0, d3.min(rows, d => d.value)) * 1.15, Math.max(0, d3.max(rows, d => d.value)) * 1.15])
    .range([F.H - M.b, M.t]);
  const x = d3.scaleBand().domain(rows.map(d => d.label))
    .range([M.l, F.W - M.r]).padding(0.3);

  forge.axLeft(svg, y, { l: M.l, n: 5, label: "count" });
  forge.axBottom(svg, x, { y: F.H - M.b, n: rows.length, tickSize: 0 });

  svg.append("g").selectAll("rect").data(rows).join("rect")
    .attr("x", d => x(d.label)).attr("width", x.bandwidth())
    .attr("y", d => Math.min(y(d.value), y(0)))
    .attr("height", d => Math.abs(y(0) - y(d.value)))
    .attr("fill", C.accent)
    .on("mousemove", (e, d) => forge.showTip(tooltip, e, `<b>${d.label}</b><br>${d.value}`, 160))
    .on("mouseout", () => forge.hideTip(tooltip));

  window.__chartReady = true;
})();
