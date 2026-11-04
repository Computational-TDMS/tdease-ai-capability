// Per-proteoform sequence coverage at backbone-bond resolution. One row per
// proteoform/hypothesis: thick band = merged cleavage-site cells (c-side /
// z-side / both), thin bands = per-spectrum tracks, ticks = modifications.
// Competing hypotheses diverge at bonds only one of them explains — the example
// ModA@K47 vs ModB@K47 demo reads differently at 14 bonds despite both being 99%.
(function () {
  const { C, tip, showTip, hideTip, axBottom } = forge;
  const D = window.DATA.data;
  const seq = D.sequence, L = seq.length, start = D.start || 1;
  const rows = D.rows;
  const CELL_NONE = 0, CELL_C = 1, CELL_Z = 2, CELL_BOTH = 3;
  const CELL_COLOR = { [CELL_C]: C.cIon, [CELL_Z]: C.zIon, [CELL_BOTH]: C.ink };
  const CELL_TEXT = { [CELL_C]: "c-side only", [CELL_Z]: "z-side only", [CELL_BOTH]: "both series" };

  const F = forge.fig({ H: 60 + rows.reduce((h, r) => h + 34 + 9 + r.tracks.length * 11 + 14, 0) + 46 });
  const { W, H, k } = F;
  const M = { l: 100 * k, r: 20 * k };
  const pw = W - M.l - M.r;
  const cw = pw / (L - 1);                 // one cell per backbone bond
  const px = r => M.l + (r - 1) * cw;      // residue r left edge (sequence numbering)
  const bondX = n => M.l + (n - 1) * cw + cw / 2;   // bond after residue n

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  /* ---------- residue ruler: position numbers every 10 + letters */
  const ruler = svg.append("g");
  ruler.append("line").attr("x1", M.l).attr("x2", W - M.r).attr("y1", 26 * k).attr("y2", 26 * k)
    .attr("stroke", C.gFaint);
  for (let p = 10 - ((start - 1) % 10); p <= L; p += 10) {
    ruler.append("text").attr("x", px(p)).attr("y", 16 * k).attr("font-size", 9 * k)
      .attr("font-family", C.fontMono).attr("fill", C.gMid).text(p + start - 1);
  }
  for (let r = 1; r <= L; r++) {
    ruler.append("text").attr("x", px(r) + cw / 2).attr("y", 40 * k).attr("font-size", 8 * k)
      .attr("font-family", C.fontMono).attr("text-anchor", "middle").attr("fill", C.gMid)
      .text(seq[r - 1])
      .append("title").text(`${seq[r - 1]}${r + start - 1}`);
  }

  /* ---------- row bands */
  let y0 = 58 * k;
  rows.forEach(row => {
    const g = svg.append("g");
    g.append("text").attr("x", M.l).attr("y", y0 + 8 * k).attr("font-size", 12 * k)
      .attr("font-weight", 700).attr("fill", C.ink).text(row.label);
    g.append("text").attr("x", W - M.r).attr("y", y0 + 8 * k).attr("font-size", 10.5 * k)
      .attr("text-anchor", "end").attr("fill", C.gMid)
      .text(`${row.coveragePct}% bonds read · ${row.nFragments} fragments`);
    y0 += 14 * k;

    /* merged band */
    drawCells(g, row.mergedSites, y0, 15 * k, row.mods);
    y0 += 17 * k;

    /* per-spectrum tracks */
    row.tracks.forEach(t => {
      g.append("text").attr("x", M.l - 8 * k).attr("y", y0 + 8 * k).attr("font-size", 9 * k)
        .attr("text-anchor", "end").attr("fill", C.gMid).text(t.name);
      drawCells(g, t.sites, y0, 7 * k, null);
      y0 += 11 * k;
    });
    y0 += 14 * k;
  });

  /* ---------- cell strip with hover + modification ticks */
  function drawCells(g, cells, y, h, mods) {
    g.append("rect").attr("x", M.l).attr("y", y).attr("width", pw).attr("height", h)
      .attr("fill", C.gGhost).attr("rx", 2);
    cells.forEach((c, i) => {
      if (!c) return;
      const r = g.append("rect")
        .attr("x", M.l + i * cw + 0.5).attr("y", y + 1 * k)
        .attr("width", Math.max(cw - 1, 1)).attr("height", h - 2 * k)
        .attr("fill", CELL_COLOR[c]).attr("rx", 1);
      r.on("mousemove", e => showTip(tt, e,
        `<b>bond after ${seq[i]}${i + start}</b><br>${CELL_TEXT[c]}`, 200))
        .on("mouseleave", () => hideTip(tt));
    });
    if (!mods) return;
    mods.forEach(m => {
      const x = m.pos === 0 ? M.l : px(m.pos) + cw / 2;
      const t = g.append("g").style("cursor", "crosshair");
      t.append("line").attr("x1", x).attr("x2", x).attr("y1", y - 12 * k).attr("y2", y + h)
        .attr("stroke", C.link).attr("stroke-width", 1.4);
      t.append("text").attr("x", x).attr("y", y - 15 * k).attr("font-size", 9 * k)
        .attr("text-anchor", "middle").attr("fill", C.link).attr("font-weight", 600)
        .text(m.pos === 0 ? `N-${m.name}` : `${seq[m.pos - 1]}${m.pos + start - 1}${m.name}`);
      t.on("mousemove", e => showTip(tt, e,
        `<b>${m.pos === 0 ? "N-terminus" : `${seq[m.pos - 1]}${m.pos + start - 1}`}</b><br>+${m.delta.toFixed(4)} Da · ${m.name}`, 210))
        .on("mouseleave", () => hideTip(tt));
    });
  }

  /* ---------- x axis + legend */
  const xa = d3.scaleLinear().domain([start, L + start - 1]).range([M.l, W - M.r]);
  axBottom(svg, xa, {
    y: H - 30 * k, n: 10,
    label: `residue position (${seq[0]}${start} … ${seq[L - 1]}${L + start - 1})`,
  });
  const legend = [["c-side (b/c ions)", C.cIon], ["z-side (y/z ions)", C.zIon],
    ["both series", C.ink], ["unread bond", C.gGhost]];
  let lx = M.l;
  legend.forEach(([label, color]) => {
    svg.append("rect").attr("x", lx).attr("y", H - 16 * k).attr("width", 10 * k).attr("height", 8 * k)
      .attr("fill", color).attr("rx", 1);
    svg.append("text").attr("x", lx + 14 * k).attr("y", H - 9 * k).attr("font-size", 10 * k)
      .attr("fill", C.gDark).text(label);
    lx += 14 * k + label.length * 5.6 * k + 22 * k;
  });
  window.__chartReady = true;
})();
