/* MS1 isotope-family map — example data (example-inputs/ms1-families.csv, the example's ladder
 * table). x = RT window centre, y = spacing (log, Da), area = member_count, colour =
 * MS2 support. The bromo ladder (Δ78.91) annotated. Family-10 sibling of feature-map. */
(function () {
  "use strict";
  const rows = window.DATA.data;
  const meta = window.DATA.meta;

  const C = forge.C;

  const data = rows.map(r => {
    const start = +r.common_rt_start, end = +r.common_rt_end;
    return {
      spacing: +r.spacing_center_da,
      rt: (start + end) / 2,
      members: +r.member_count,
      intensity: +r.intensity_total,
      ms2: r.ms2_support && r.ms2_support !== "none",
      masses: r.masses.split(";").map(s => +s),
    };
  });

  const W = 1054, H = 620, M = { t: 20, r: 30, b: 52, l: 70 };
  const svg = d3.select("#chart").append("svg")
    .attr("viewBox", `0 0 ${W} ${H}`).attr("width", W).attr("height", H)
    .style("font-family", C.font);

  const x = d3.scaleLinear().domain(d3.extent(data, d => d.rt)).nice().range([M.l, W - M.r]);
  const y = d3.scaleLog().domain(d3.extent(data, d => d.spacing)).nice().range([H - M.b, M.t]);
  const rSize = d3.scaleSqrt().domain(d3.extent(data, d => d.members)).range([3, 13]);

  x.ticks(7).forEach(v => svg.append("line").attr("x1", x(v)).attr("x2", x(v))
    .attr("y1", M.t).attr("y2", H - M.b).attr("stroke", C.grid).attr("stroke-width", 0.5));
  y.ticks(5).forEach(v => svg.append("line").attr("x1", M.l).attr("x2", W - M.r)
    .attr("y1", y(v)).attr("y2", y(v)).attr("stroke", C.grid).attr("stroke-width", 0.5));

  const tooltip = forge.tip();

  svg.append("g").selectAll("circle").data(data).join("circle")
    .attr("cx", d => x(d.rt)).attr("cy", d => y(d.spacing))
    .attr("r", d => rSize(d.members))
    .attr("fill", d => (d.ms2 ? C.accent : C.unmod)).attr("opacity", 0.7)
    .attr("stroke", "#fff").attr("stroke-width", 0.8)
    .on("mousemove", (e, d) => forge.showTip(tooltip, e,
      `<b>Δ${d.spacing.toFixed(2)} Da</b> family · ${d.members} members · RT ${d.rt.toFixed(1)} min<br>` +
      `masses: ${d.masses.slice(0, 4).map(m => m.toFixed(1)).join(" → ")}${d.masses.length > 4 ? " → …" : ""}<br>` +
      `<span style="opacity:.75">MS2 support: ${d.ms2 ? "yes" : "none"} · total intensity ${d.intensity.toExponential(1)}</span>`, 460))
    .on("mouseout", () => forge.hideTip(tooltip));

  /* annotate the bromo ladder family (the example's story lives here) */
  const bromo = data.find(d => Math.abs(d.spacing - 78.91) < 0.3);
  if (bromo) {
    svg.append("circle").attr("cx", x(bromo.rt)).attr("cy", y(bromo.spacing))
      .attr("r", rSize(bromo.members) + 4).attr("fill", "none")
      .attr("stroke", C.accent).attr("stroke-width", 1.2).attr("stroke-dasharray", "3 2");
    svg.append("text").attr("x", x(bromo.rt)).attr("y", y(bromo.spacing) - rSize(bromo.members) - 8)
      .attr("text-anchor", "middle").attr("font-size", 10.5).attr("fill", C.ink)
      .attr("paint-order", "stroke").attr("stroke", "#fff").attr("stroke-width", 2.5)
      .text("bromo 78.91");
  }

  svg.append("g").attr("transform", `translate(0,${H - M.b})`)
    .call(d3.axisBottom(x).ticks(7).tickSize(3))
    .call(s => s.select(".domain").attr("stroke-width", 0.75));
  svg.append("g").attr("transform", `translate(${M.l},0)`)
    .call(d3.axisLeft(y).ticks(5).tickSize(3))
    .call(s => s.select(".domain").attr("stroke-width", 0.75));
  svg.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", H - 10)
    .attr("text-anchor", "middle").attr("font-size", 13).text("family RT window centre (min)");
  svg.append("text").attr("transform", "rotate(-90)").attr("x", -(M.t + H - M.b) / 2)
    .attr("y", 16).attr("text-anchor", "middle").attr("font-size", 13)
    .text("isotope spacing (Da, log)");

  forge.legend(svg, [
    { label: "MS2-supported family", color: C.accent },
    { label: "no MS2 support yet", color: C.unmod },
  ], { x: M.l + 12, y: M.t + 13, pitch: 150, fs: 10.5 });
  svg.append("text").attr("x", M.l + 310).attr("y", M.t + 16.5).attr("font-size", 10.5)
    .attr("fill", C.gMid).text("area = member count");

  window.__chartReady = true;
})();
