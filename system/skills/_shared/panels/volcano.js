;(function () {
  "use strict";
  const data = window.DATA.data || {};
  const points = data.points || [];
  const C = forge.C;
  const { W, H } = forge.fig({ W: 900, H: 640 });
  const M = { t: 48, r: 190, b: 70, l: 82 };
  const plotW = W - M.l - M.r, plotH = H - M.t - M.b;
  const svg = d3.select("#chart").append("svg")
    .attr("viewBox", `0 0 ${W} ${H}`).attr("width", W).attr("height", H)
    .style("font-family", C.font);
  const extent = values => {
    const e = d3.extent(values.filter(Number.isFinite));
    return e[0] === e[1] ? [e[0] - 1, e[1] + 1] : e;
  };
  const fcCutoff = Number.isFinite(+data.fcCutoff) ? +data.fcCutoff : 1;
  const fdr = Number.isFinite(+data.fdr) ? +data.fdr : 0.05;
  const thresholdY = -Math.log10(fdr);
  const xMax = Math.max(fcCutoff * 1.25, d3.max(points, d => Math.abs(+d.x)) || 1) * 1.05;
  const yMax = Math.max(thresholdY * 1.15, d3.max(points, d => +d.y) || 1) + 0.35;
  const x = d3.scaleLinear().domain([-xMax, xMax]).range([M.l, M.l + plotW]);
  const y = d3.scaleLinear().domain([0, yMax]).range([M.t + plotH, M.t]);
  const sizeValues = points.map(d => +d.sizeValue).filter(Number.isFinite);
  const sizeExtent = d3.extent(sizeValues);
  const sizeScale = sizeValues.length && sizeExtent[0] !== sizeExtent[1]
    ? d3.scaleSqrt().domain(sizeExtent).range([16, 78]) : null;
  const radius = d => Math.sqrt((sizeScale ? sizeScale(+d.sizeValue) : 24) / Math.PI);
  const color = {
    ns: C.ns, "fc-only": C.fcOnly, "p-only": C.pOnly,
    significant: C.accent
  };
  const label = {
    ns: "not significant", "fc-only": "|log2FC| only", "p-only": "BH FDR only",
    significant: "BH FDR + |log2FC|"
  };
  const esc = value => String(value ?? "").replace(/[&<>\"]/g, c => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;"
  }[c]));
  const tip = forge.tip();
  forge.gridV(svg, x, { h: plotH, t: M.t, n: 7 });
  forge.gridH(svg, y, { w: plotW, l: M.l, n: 6 });
  const tint = (x0, x1) => svg.append("rect")
    .attr("x", x(x0)).attr("y", M.t).attr("width", x(x1) - x(x0))
    .attr("height", Math.max(0, y(thresholdY) - M.t))
    .attr("fill", C.accent).attr("fill-opacity", 0.045);
  tint(fcCutoff, x.domain()[1]);
  tint(x.domain()[0], -fcCutoff);
  [-fcCutoff, fcCutoff].forEach(v => svg.append("line")
    .attr("x1", x(v)).attr("x2", x(v)).attr("y1", M.t).attr("y2", M.t + plotH)
    .attr("stroke", C.threshold).attr("stroke-dasharray", "6 4"));
  svg.append("line").attr("x1", M.l).attr("x2", M.l + plotW)
    .attr("y1", y(thresholdY)).attr("y2", y(thresholdY))
    .attr("stroke", C.threshold).attr("stroke-dasharray", "6 4");
  svg.append("text").attr("x", M.l + plotW - 4).attr("y", y(thresholdY) - 7)
    .attr("text-anchor", "end").attr("font-size", 9).attr("fill", C.gMid)
    .text(`BH FDR ${d3.format(".2g")(fdr)} · |FC| ${d3.format(".2g")(fcCutoff)}`);

  const marks = svg.append("g").selectAll("circle").data(points).join("circle")
    .attr("cx", d => x(+d.x)).attr("cy", d => y(+d.y)).attr("r", radius)
    .attr("fill", d => color[d.state] || C.ns)
    .attr("fill-opacity", d => d.state === "ns" ? 0.42 : 0.86)
    .attr("stroke", "#fff").attr("stroke-width", 0.7).style("cursor", "pointer")
    .on("mouseenter", function (event, d) {
      d3.select(this).attr("stroke", C.ink).attr("stroke-width", 1.8);
      const capped = d.capped ? "<br><em>y capped at display limit</em>" : "";
      const html = `<b>${esc(d.feature)}</b><br>log2FC ${d3.format("+.3f")(+d.x)}` +
        `<br>BH FDR ${d.fdr == null ? "—" : d3.format(".3e")(+d.fdr)}` +
        `<br>${esc(label[d.state] || d.state)} · mean abundance ${d.sizeValue == null ? "—" : d3.format(",.1f")(+d.sizeValue)}` + capped;
      forge.showTip(tip, event, html, 310);
    })
    .on("mousemove", (event) => forge.showTip(tip, event, tip.html(), 310))
    .on("mouseleave", function (_, d) {
      d3.select(this).attr("stroke", "#fff").attr("stroke-width", 0.7);
      forge.hideTip(tip);
    })
    .on("click", (_, d) => window.reportSelect("volcano", {
      kind: "feature", featureId: d.feature, log2fc: d.x, fdr: d.fdr
    }));

  const labeled = points.filter(d => d.state === "significant")
    .sort((a, b) => (+a.fdr || 1) - (+b.fdr || 1)).slice(0, 4);
  const levels = forge.dodgeX(labeled.map(d => x(+d.x)), 72, 3);
  const labelLayer = svg.append("g");
  labeled.forEach((d, i) => {
    const px = x(+d.x), py = y(+d.y), side = px < M.l + plotW / 2 ? -1 : 1;
    const level = levels[i] < 0 ? 0 : levels[i];
    const tx = px + side * (44 + level * 18), ty = Math.max(M.t + 12, py - 18 - level * 18);
    labelLayer.append("path").attr("d", `M${px + side * 4},${py - 2} L${tx - side * 4},${ty + 3}`)
      .attr("fill", "none").attr("stroke", C.accent).attr("stroke-width", 0.8);
    labelLayer.append("text").attr("x", tx).attr("y", ty).attr("text-anchor", side < 0 ? "end" : "start")
      .attr("font-size", 9).attr("font-weight", 700).attr("fill", C.accent)
      .text(String(d.feature).slice(0, 24));
  });
  if (labeled.length) svg.append("text").attr("x", M.l + 8).attr("y", M.t + 14)
    .attr("font-size", 9).attr("fill", C.gMid).text(`labels: top ${labeled.length} thresholded features`);

  forge.axBottom(svg, x, { y: M.t + plotH, n: 7,
    label: data.xlabel || "log2 fold change" });
  forge.axLeft(svg, y, { l: M.l, n: 6,
    label: data.ylabel || "-log10 BH FDR" });

  const legend = svg.append("g").attr("transform", `translate(${M.l + plotW + 22},${M.t})`);
  legend.append("text").attr("font-size", 11).attr("font-weight", 700)
    .attr("fill", C.ink).text("Point state");
  const counts = d3.rollup(points, v => v.length, d => d.state);
  ["significant", "p-only", "fc-only", "ns"].forEach((state, i) => {
    const row = legend.append("g").attr("transform", `translate(0,${22 + i * 25})`);
    row.append("circle").attr("cx", 5).attr("cy", 0).attr("r", 5)
      .attr("fill", color[state]).attr("fill-opacity", state === "ns" ? 0.42 : 0.86);
    row.append("text").attr("x", 17).attr("y", 4).attr("font-size", 9)
      .attr("fill", C.gMid).text(`${label[state]} (${counts.get(state) || 0})`);
  });
  if (sizeScale) {
    const sy = 36 + 4 * 25;
    legend.append("text").attr("y", sy).attr("font-size", 11).attr("font-weight", 700)
      .attr("fill", C.ink).text(data.sizeLabel || "mean abundance");
    [sizeExtent[0], d3.median(sizeValues), sizeExtent[1]].forEach((value, i) => {
      const r = radius({ sizeValue: value });
      const row = legend.append("g").attr("transform", `translate(${r + 4},${sy + 28 + i * 28})`);
      row.append("circle").attr("cy", -r).attr("r", r).attr("fill", C.gSoft)
        .attr("fill-opacity", 0.22).attr("stroke", C.gMid).attr("stroke-width", 0.7);
      row.append("text").attr("x", 18).attr("y", 4 - r).attr("font-size", 9)
        .attr("fill", C.gMid).text(d3.format(",.1f")(value));
    });
  }
  forge.hideTip(tip);
  window.__chartReady = true;
})();
