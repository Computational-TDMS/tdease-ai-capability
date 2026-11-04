// Report panel: vertical bars, optionally grouped by series.
// Contract: window.DATA.data = { rows: [{label, value, err?, points?, hover?, series?}],
//   series?: ["A","B"], xlabel?, ylabel?, valueLabels?: bool, height?: px,
//   colors?: "by-row" | "vik" | "batlow" }
// `err` (same units as value) draws a ± whisker around the bar top;
// `points` (individual observations) are drawn as dodged dots over the bar;
// `colors` recolors ungrouped rows: "by-row" cycles the categorical palette,
// a ramp name samples that crameri ramp across the rows (vik: blue↔pink,
// skipping its near-white middle; a single row takes the ramp's mid-blue).
// Bars start at zero (counts/percentages); the axis is a view, the CSV/report
// envelope is the record.
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, gridH } = forge;
  const { rows, series, xlabel, ylabel, valueLabels, height, colors } = window.DATA.data;
  const F = forge.fig({ W: 1054, H: height || 300 });
  const { W, H, k } = F;
  const M = { l: 58 * k, r: 22 * k, t: 26 * k, b: 62 * k };

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip();

  const labels = [...new Set(rows.map(r => r.label))];
  const sers = series || [null];
  const ramp = colors === "vik" && typeof vik !== "undefined" ? vik
    : colors === "batlow" && typeof batlow !== "undefined" ? batlow : null;
  // Sample a ramp across k rows: blue arm then red/pink arm, skipping the
  // diverging near-white middle and vik's terracotta band.
  const rampColor = (k, i) => {
    if (k === 1) return ramp[Math.round(0.25 * (ramp.length - 1))];
    const t = i / (k - 1);
    const pos = t <= 0.5 ? 0.05 + (t / 0.5) * 0.35 : 0.62 + ((t - 0.5) / 0.5) * 0.34;
    return ramp[Math.round(pos * (ramp.length - 1))];
  };
  const colorOf = r => r.series
    ? C.cat[sers.indexOf(r.series) % C.cat.length]
    : ramp ? rampColor(rows.length, rows.indexOf(r))
    : colors === "by-row" ? C.cat[rows.indexOf(r) % C.cat.length]
    : C.cat[0];
  const x = d3.scaleBand().domain(labels).range([M.l, W - M.r]).paddingInner(0.25);
  const x1 = sers.length > 1
    ? d3.scaleBand().domain(sers).range([0, x.bandwidth()]).paddingInner(0.15)
    : null;
  const yMax = Math.max(1e-9, d3.max(rows, r =>
    Math.max(r.value + (r.err || 0), ...(r.points || [0]))));
  const y = d3.scaleLinear().domain([0, yMax * 1.12]).range([H - M.b, M.t]);
  const band = x1 ? x1.bandwidth() : x.bandwidth();
  const bx = r => x(r.label) + (x1 ? x1(r.series || sers[0]) : 0);

  gridH(svg.append("g"), y, { w: W - M.l - M.r, l: M.l, n: 5 });

  svg.append("g").selectAll("rect").data(rows).join("rect")
    .attr("x", bx)
    .attr("width", band)
    .attr("y", r => y(r.value))
    .attr("height", r => y(0) - y(r.value))
    .attr("fill", colorOf)
    .style("cursor", "crosshair")
    .on("mousemove", (e, r) => showTip(tt, e, r.hover ||
      `<b>${r.label}</b>${r.series ? " · " + r.series : ""}: ${d3.format(",.6~g")(r.value)}`))
    .on("mouseleave", hideTip(tt));

  // ± whiskers for rows carrying err; clamped at zero.
  const errs = rows.filter(r => r.err > 0);
  if (errs.length) {
    const cx = r => bx(r) + band / 2;
    const cap = Math.min(9 * k, band / 3);
    svg.append("g").selectAll("line").data(errs).join("line")
      .attr("x1", cx).attr("x2", cx)
      .attr("y1", r => y(Math.max(0, r.value - r.err))).attr("y2", r => y(r.value + r.err))
      .attr("stroke", C.gDark).attr("stroke-width", 1);
    for (const yy of [r => y(r.value + r.err), r => y(Math.max(0, r.value - r.err))]) {
      svg.append("g").selectAll("line").data(errs).join("line")
        .attr("x1", r => cx(r) - cap).attr("x2", r => cx(r) + cap)
        .attr("y1", yy).attr("y2", yy)
        .attr("stroke", C.gDark).attr("stroke-width", 1);
    }
  }

  // Individual observations as dodged dots over the bar.
  const pts = rows.flatMap(r => (r.points || []).map((v, i) => ({ r, v, i })));
  if (pts.length) {
    const span = Math.min(16 * k, band / 2.2);
    svg.append("g").selectAll("circle").data(pts).join("circle")
      .attr("cx", d => {
        const n = d.r.points.length;
        const step = n > 1 ? span / (n - 1) : 0;
        return bx(d.r) + band / 2 + (d.i - (n - 1) / 2) * step;
      })
      .attr("cy", d => y(Math.max(0, d.v)))
      .attr("r", 2.7 * k)
      .attr("fill", C.gDark)
      .style("cursor", "crosshair")
      .on("mousemove", (e, d) => showTip(tt, e,
        `<b>${d.r.label}</b>: ${d3.format(",.6~g")(d.v)}`))
      .on("mouseleave", hideTip(tt));
  }

  if (valueLabels && (!x1 || labels.length <= 16)) {
    svg.append("g").selectAll("text").data(rows).join("text")
      .attr("x", r => bx(r) + band / 2)
      .attr("y", r => y(r.value + (r.err || 0)) - 4)
      .attr("text-anchor", "middle").attr("font-size", 9.5 * k)
      .attr("fill", C.gDark)
      .text(r => d3.format(",.6~g")(r.value));
  }

  axBottom(svg.append("g"), x, { y: y(0), n: Math.min(labels.length, 12), label: xlabel });
  axLeft(svg.append("g"), y, { l: M.l, n: 5, label: ylabel });

  // Long or numerous labels rotate; short few stay horizontal.
  if (labels.some(l => String(l).length > 8) || labels.length > 10) {
    svg.selectAll(".tick text").attr("transform", "rotate(-45)")
      .attr("text-anchor", "end").attr("dx", "-0.4em").attr("dy", "0.65em");
  }

  if (sers.length > 1) {
    const legend = svg.append("g").attr("transform", `translate(${W - M.r - 10 * sers.length * 7 * k},${M.t - 12 * k})`);
    sers.forEach((s, i) => {
      const g = legend.append("g").attr("transform", `translate(${i * 84 * k},0)`);
      g.append("rect").attr("width", 10).attr("height", 10)
        .attr("fill", C.cat[sers.indexOf(s) % C.cat.length]).attr("rx", 2);
      g.append("text").attr("x", 14).attr("y", 9).attr("font-size", 10.5 * k)
        .attr("fill", C.gDark).text(s);
    });
  }
})();
