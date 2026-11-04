// Report panel: heatmap matrix (sample correlation, expression heatmap).
// Contract: window.DATA.data = {
//   rowLabels: [..], colLabels: [..], values: [[number|null]],
//   diverging?: bool, lo: number, hi: number,     // explicit — deterministic, computed python-side
//   colGroups?: [{label, n}],                     // vertical separators between group blocks
//   legendTitle?: string, xlabel?, ylabel?, height?: px, fmt?: ".2f" }
// null = structural absence (not detected) — distinct grey, never the low colour.
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, cb } = forge;
  const { rowLabels, colLabels, values, diverging, lo, hi, colGroups,
          legendTitle, xlabel, ylabel, height, fmt = ".2f" } = window.DATA.data;
  const nR = rowLabels.length, nC = colLabels.length;
  const F = forge.fig({ W: 1054, H: height || Math.max(340, 90 + nR * 16) });
  const { W, H, k } = F;
  const M = { l: 132 * k, r: 26 * k, t: 30 * k, b: 96 * k };

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip();

  const ramp = diverging ? vik : batlow;
  // piecewise through the full Crameri ramp; diverging ramps are symmetric,
  // so norm(lo..hi) over [0,1] puts the pale midpoint at the domain centre.
  const norm = d3.scaleLinear().domain([lo, hi]).range([0, 1]).clamp(true);
  const piecewise = d3.piecewise(d3.interpolateRgb, ramp);
  const color = v => piecewise(norm(v));

  const w = W - M.l - M.r, h = H - M.t - M.b;
  const cw = w / nC, ch = h / nR;
  const xPos = j => M.l + j * cw;
  const yPos = i => M.t + i * ch;
  const format = d3.format(fmt);

  const cells = [];
  values.forEach((row, i) => row.forEach((v, j) => cells.push({ i, j, v })));

  svg.append("g").selectAll("rect").data(cells).join("rect")
    .attr("x", c => xPos(c.j)).attr("y", c => yPos(c.i))
    .attr("width", cw).attr("height", ch)
    .attr("fill", c => c.v === null ? C.gFaint : color(c.v))
    .on("mousemove", (e, c) => showTip(tt, e,
      `<b>${rowLabels[c.i]}</b> × <b>${colLabels[c.j]}</b><br>` +
      (c.v === null ? "not detected" : format(c.v))))
    .on("mouseleave", () => hideTip(tt));

  // group block separators (pre-ordered python-side; lines make the design readable)
  if (colGroups) {
    let acc = 0;
    for (const g of colGroups.slice(0, -1)) {
      acc += g.n;
      svg.append("line").attr("x1", xPos(acc)).attr("x2", xPos(acc))
        .attr("y1", M.t).attr("y2", M.t + h)
        .attr("stroke", C.ink).attr("stroke-width", 1.4);
    }
  }

  // scaleBand, not scalePoint: d3-axis emits no ticks for point scales (v7).
  axBottom(svg.append("g"), d3.scaleBand(colLabels, [M.l, M.l + w]),
    { y: M.t + h, n: nC, label: xlabel });
  svg.selectAll(".tick text").attr("transform", "rotate(-50)")
    .attr("text-anchor", "end").attr("dx", "-0.3em").attr("dy", "0.55em");
  axLeft(svg.append("g"), d3.scaleBand(rowLabels, [M.t, M.t + h]),
    { l: M.l, n: Math.min(nR, 20), label: ylabel });
  // Row labels must fit inside a cell: cap at ~55% row height, k-normalized.
  // Applies to both axes' ticks (the rotated column labels are smaller anyway).
  svg.selectAll(".tick text").attr("font-size", Math.min(10, (ch * 0.55) / k + 5));

  cb(svg, color, { lo, hi, W, y: H - 30, title: legendTitle,
                   fmt: format });
})();
