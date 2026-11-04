// Proteoform map — every proteoform of one accession as a packed horizontal
// sequence bar (non-overlapping intervals share a lane, short forms tile).
// Colour = PrSM count (batlow reversed, darker = more). Identified
// modification families ride ON the bars as numbered boxes. The figure legend
// is ONE right-hand column, matching the plotly original's layout: numbered
// family notes with leader lines, then the vertical PrSM colourbar, then the
// mark legend. Lane packing and mod placement are computed in gen-data.py
// (ports of tdp-opensearch figures.py).
(function () {
  const { C, tip, showTip, hideTip, axBottom } = forge;
  const D = window.DATA.data;
  const rows = D.rows, L = D.proteinLength, nLanes = D.nLanes, boxes = D.boxes;

  const H_DESIGN = Math.max(
    34 + nLanes * 13 + 64,                        // plot: lanes + parent strip + axis
    34 + boxes.length * 58 + 18 + 104 + 28 + 44,  // right column: notes + colourbar + marks
  );
  const F = forge.fig({ H: H_DESIGN });
  const { W, H, k } = F;
  const M = { l: 64 * k };
  const NOTE_W = 264 * k;
  const PLOT_R = W - NOTE_W - 24 * k;
  const RX = PLOT_R + 28 * k;                     // right legend column x
  const top = 34 * k;
  const laneY = l => top + l * 13 * k;
  const barH = 7 * k;
  const yParent = laneY(nLanes) + 10 * k;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  const x = d3.scaleLinear().domain([0.5, L + 0.5]).range([M.l, PLOT_R]);
  const maxCount = d3.max(rows, r => r.count) || 1;
  const ramp = d3.scaleSequential(d3.interpolateRgbBasis([...batlow].reverse()))
    .domain([0, Math.sqrt(maxCount)]);
  const barColor = c => ramp(Math.sqrt(c));

  /* ---------- bars + mod marks */
  rows.forEach(r => {
    const yC = laneY(r.lane) + barH / 2;
    const g = svg.append("g");
    g.append("rect")
      .attr("x", x(r.first)).attr("y", laneY(r.lane))
      .attr("width", Math.max(x(r.last) - x(r.first), 1.5)).attr("height", barH)
      .attr("rx", 1.5).attr("fill", barColor(r.count))
      .on("mousemove", e => showTip(tt, e,
        `<b>${D.accession} · residues ${r.first}\u2013${r.last}</b><br>` +
        `${r.count} PrSM${r.count === 1 ? "" : "s"} · lane ${r.lane + 1}<br>` +
        `<span style="font-family:${C.fontMono};font-size:10px">${r.proforma.slice(0, 80)}${r.proforma.length > 80 ? "…" : ""}</span>`,
        340))
      .on("mouseleave", () => hideTip(tt));
    r.mods.forEach(m => {
      const mx = x(r.first + m.pos - 1);
      if (m.box) {
        const note = boxes.find(b => b.n === m.box);
        const bg = svg.append("g").style("cursor", "crosshair");
        bg.append("rect")
          .attr("x", mx - 5.5 * k).attr("y", yC - 5.5 * k)
          .attr("width", 11 * k).attr("height", 11 * k)
          .attr("fill", "#fff").attr("stroke", C.ink).attr("stroke-width", 1.2);
        bg.append("text").attr("x", mx).attr("y", yC + 3 * k)
          .attr("text-anchor", "middle").attr("font-size", 8.5 * k)
          .attr("font-weight", 700).attr("fill", C.ink).text(m.box);
        bg.on("mousemove", e => showTip(tt, e,
          `<b>#${m.box} ${note.label}</b><br>+${m.mass.toFixed(4)} Da · on ${note.nBars} bars`, 300))
          .on("mouseleave", () => hideTip(tt));
      } else {
        svg.append("circle").attr("cx", mx).attr("cy", yC).attr("r", 2.4 * k)
          .attr("fill", C.gMid)
          .on("mousemove", e => showTip(tt, e,
            m.name ? `<b>${m.name}</b> @${m.pos}` : `<b>+${m.mass.toFixed(4)} Da</b> · unresolved`, 200))
          .on("mouseleave", () => hideTip(tt));
      }
    });
  });

  /* ---------- parent sequence strip + axis */
  svg.append("rect").attr("x", x(1)).attr("y", yParent).attr("width", x(L) - x(1))
    .attr("height", 8 * k).attr("rx", 2).attr("fill", C.gSoft);
  svg.append("text").attr("x", x(1)).attr("y", yParent - 5 * k)
    .attr("font-size", 10 * k).attr("fill", C.gMid).text(`parent sequence · ${L} aa`);
  axBottom(svg, x, { y: yParent + 10 * k, n: 10, label: `residue position (${D.accession})` });

  /* ---------- right legend column ①: numbered family notes + leader lines */
  const noteX = RX;
  boxes.forEach((b, i) => {
    const ny = top + 6 * k + i * 58 * k;
    const strongest = rows.filter(r => r.mods.some(m => m.box === b.n))
      .reduce((a, r) => (r.count > a.count ? r : a), rows[0]);
    const mod = strongest.mods.find(m => m.box === b.n);
    const tx = x(strongest.first + mod.pos - 1), ty = laneY(strongest.lane) + barH / 2;
    svg.append("path")
      .attr("d", `M${tx + 6 * k},${ty} L${tx + 14 * k},${ty} L${noteX - 10 * k},${ny + 6 * k}`)
      .attr("fill", "none").attr("stroke", C.gGhost).attr("stroke-width", 1.1);
    svg.append("rect").attr("x", noteX).attr("y", ny).attr("width", 11 * k).attr("height", 11 * k)
      .attr("fill", "#fff").attr("stroke", C.ink).attr("stroke-width", 1.2);
    svg.append("text").attr("x", noteX + 5.5 * k).attr("y", ny + 8.5 * k)
      .attr("text-anchor", "middle").attr("font-size", 8.5 * k)
      .attr("font-weight", 700).attr("fill", C.ink).text(b.n);
    const words = b.label.split(" ");
    const lines = [];
    for (let w = 0; w < words.length;) {
      let line = "";
      while (w < words.length && (line + " " + words[w]).trim().length <= 30) line = (line + " " + words[w++]).trim();
      lines.push(line);
    }
    lines.slice(0, 2).forEach((ln, li) => svg.append("text")
      .attr("x", noteX + 17 * k).attr("y", ny + 8 * k + li * 12 * k)
      .attr("font-size", 10 * k).attr("font-weight", li === 0 ? 600 : 400)
      .attr("fill", li === 0 ? C.ink : C.gDark).text(ln + (li === 1 && lines.length > 2 ? "…" : "")));
    const masses = b.masses.length > 3
      ? `+${b.masses[0].toFixed(2)} … +${b.masses[b.masses.length - 1].toFixed(2)} Da`
      : b.masses.map(m => `+${m.toFixed(2)}`).join(" / ") + " Da";
    svg.append("text").attr("x", noteX + 17 * k).attr("y", ny + 8 * k + Math.min(lines.length, 2) * 12 * k)
      .attr("font-size", 9 * k).attr("font-family", C.fontMono).attr("fill", C.gMid)
      .text(`${masses} · ${b.nBars} bars`);
  });

  /* ---------- right legend column ②: vertical PrSM colourbar */
  const cbY = top + 6 * k + Math.max(boxes.length, 2) * 58 * k;
  const cbW = 13 * k, cbH = 104 * k;
  const defs = svg.append("defs");
  const grad = defs.append("linearGradient").attr("id", "pf-cb")
    .attr("x1", 0).attr("y1", 1).attr("x2", 0).attr("y2", 0);
  d3.range(21).forEach(i => grad.append("stop")
    .attr("offset", i * 5 + "%")
    .attr("stop-color", barColor(maxCount * i / 20)));
  svg.append("text").attr("x", RX).attr("y", cbY - 7 * k)
    .attr("font-size", 10 * k).attr("font-weight", 600).attr("fill", C.gDark)
    .text("PrSM count");
  svg.append("rect").attr("x", RX).attr("y", cbY).attr("width", cbW).attr("height", cbH)
    .attr("fill", "url(#pf-cb)").attr("rx", 2);
  svg.append("text").attr("x", RX + cbW + 7 * k).attr("y", cbY + 8 * k)
    .attr("font-size", 9 * k).attr("font-family", C.fontMono).attr("fill", C.gMid)
    .text(maxCount.toFixed(0));
  svg.append("text").attr("x", RX + cbW + 7 * k).attr("y", cbY + cbH)
    .attr("font-size", 9 * k).attr("font-family", C.fontMono).attr("fill", C.gMid)
    .text("0");

  /* ---------- right legend column ③: mark legend */
  let ly = cbY + cbH + 28 * k;
  svg.append("circle").attr("cx", RX + 4 * k).attr("cy", ly - 3 * k).attr("r", 3 * k)
    .attr("fill", C.gMid);
  svg.append("text").attr("x", RX + 13 * k).attr("y", ly)
    .attr("font-size", 10 * k).attr("fill", C.gDark).text("unresolved shift");
  ly += 24 * k;
  svg.append("rect").attr("x", RX - 1 * k).attr("y", ly - 9 * k)
    .attr("width", 11 * k).attr("height", 11 * k)
    .attr("fill", "#fff").attr("stroke", C.ink).attr("stroke-width", 1.2);
  svg.append("text").attr("x", RX + 4 * k).attr("y", ly - 1 * k)
    .attr("text-anchor", "middle").attr("font-size", 8.5 * k)
    .attr("font-weight", 700).attr("fill", C.ink).text("n");
  svg.append("text").attr("x", RX + 13 * k).attr("y", ly)
    .attr("font-size", 10 * k).attr("fill", C.gDark)
    .text("identified family → note above");
  window.__chartReady = true;
})();
