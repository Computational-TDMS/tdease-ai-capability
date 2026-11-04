/* chart-forge component: MS/MS annotation panels.
 * Distilled from the frontend AnnotatedSpectrumIPSA.vue (IPSA / TopPIC-PRSM /
 * TopMSV lineage) into plain D3 + forge tokens. Inlined by build.py --script
 * BEFORE the chart script; requires window.forge (templates/micro.js).
 *
 * Exposes window.forgeMs2:
 *   T            PRSM title geometry constants (30 AA rows, 10-AA blocks)
 *   ionColor(ion)         a/b/c/x/y/z -> token colour
 *   fragmentLabel(f)      "c18^2+-H2O" convention
 *   fkey(f)               group key ion+number+charge
 *   cleavageOf(f, L)      residue index the cleavage sits AFTER — topological
 *                         truth: c_n cleaves after n; z_n after L-n. (The
 *                         frontend drew suffix ions at `position` directly,
 *                         which only works for pre-converted data.)
 *   prsmTitle(svg, residues, frags, o) -> {height, xOf, yOf, marks, letters}
 *   errorPanel(svg, frags, o)          -> {dots, yOf}
 *   annBands(g, anns, o)     standard `annotations` contract: band layer
 *   annCallouts(g, anns, o)  standard `annotations` contract: callout layer
 *
 * Consumed by examples/fragment-annotation, examples/mgf-accumulated;
 * add consumers before growing API.
 */
window.forgeMs2 = (function () {
  "use strict";

  const T = {
    rowLen: 30, blockLen: 10,          // TopPIC PRSM: 30 residues/row, blocks of 10
    letterW: 28, gapW: 20, rowH: 40,
    top: 46, left: 50, numW: 20,
    letterFs: 16,
  };

  function ionColor(ion) {
    const C = forge.C;
    const m = { a: C.aIon, b: C.bIon, c: C.cIon, x: C.xIon, y: C.yIon, z: C.zIon };
    return m[ion] || C.gMid;
  }

  const fragmentLabel = f =>
    f.ion + f.number + (f.charge > 1 ? "^" + f.charge + "+" : "")
      + (f.loss ? "-" + f.loss : "");

  const fkey = f => f.ion + f.number + "+" + f.charge + (f.loss ? "-" + f.loss : "");

  const cleavageOf = (f, L) => (f.ion === "a" || f.ion === "b" || f.ion === "c")
    ? f.number : L - f.number;

  /* modification tag -> display text: numeric masses get an explicit sign */
  function modText(m) {
    const n = parseFloat(m);
    return Number.isFinite(n) ? (n >= 0 ? "+" : "") + n.toFixed(2) : String(m);
  }

  /* ---- PRSM title: multi-row annotated sequence with cleavage polylines ----
   * residues: [{pos, char, mod}] (mod = tag string or null)
   * frags:    matched fragments; one polyline per ion+number group (strongest
   *            assignment wins), suffix ions mirrored below the row.
   * o: {top, onResidue(e, pos), onResidueLeave, onCleavage(e, frag), onCleavageLeave}
   * Craft rules (frontend-faithful): cleavage marks never cross a row break
   * and never sit adjacent to a modified residue — both would overprint. */
  function prsmTitle(svg, residues, frags, o) {
    const C = forge.C;
    const L = residues.length;
    const byPos = new Map(residues.map(r => [r.pos, r]));
    const k = o.k || 1;                     // journal builds shrink the type grid
    const rowLen = T.rowLen, blockLen = T.blockLen;
    const letterW = T.letterW * k, gapW = T.gapW * k, rowH = T.rowH * k;
    const top = T.top * k, left = T.left * k, numW = T.numW * k;
    const letterFs = Math.max(T.letterFs * k, 9);

    const rows = [];
    for (let s = 1; s <= L; s += rowLen) {
      rows.push({ start: s, end: Math.min(s + rowLen - 1, L) });
    }
    const rowOf = pos => rows.findIndex(r => pos >= r.start && pos <= r.end);
    const xOf = pos => {
      const inRow = pos - rows[rowOf(pos)].start;
      return left + numW + inRow * letterW + Math.floor(inRow / blockLen) * gapW;
    };
    const yOf = pos => top + rowOf(pos) * rowH;
    const height = top + (rows.length - 1) * rowH + 34 * k;

    const g = svg.append("g").attr("transform", `translate(0,${o.top || 0})`)
      .attr("font-family", C.fontMono);

    /* position numbers at each row's start and end */
    rows.forEach(r => {
      if (r.end < r.start) return;
      g.append("text").attr("x", left).attr("y", yOf(r.start))
        .attr("text-anchor", "end").attr("font-size", Math.max(12 * k, 9))
        .attr("fill", C.gMid).text(r.start);
      g.append("text").attr("x", xOf(r.end) + numW + 12 * k).attr("y", yOf(r.end))
        .attr("font-size", Math.max(12 * k, 9)).attr("fill", C.gMid).text(r.end);
    });

    /* residues: modified sites get a background marker + mass label above */
    const letters = g.selectAll(null).data(residues).join("text")
      .attr("x", d => xOf(d.pos)).attr("y", d => yOf(d.pos))
      .attr("text-anchor", "middle").attr("dominant-baseline", "middle")
      .attr("font-size", letterFs)
      .attr("font-weight", d => (d.mod ? 700 : 500))
      .attr("fill", d => (d.mod ? C.accent : C.ink))
      .style("cursor", "default")
      .text(d => d.char)
      .on("mouseenter", (e, d) => o.onResidue && o.onResidue(e, d.pos))
      .on("mouseleave", () => o.onResidueLeave && o.onResidueLeave());

    residues.forEach(d => {
      if (!d.mod) return;
      g.append("rect")
        .attr("x", xOf(d.pos) - letterW * 0.45).attr("y", yOf(d.pos) - letterFs * 0.75)
        .attr("width", letterW * 0.9).attr("height", letterFs * 1.4).attr("rx", 2)
        .attr("fill", C.accent).attr("fill-opacity", 0.12);
      g.append("text")
        .attr("x", xOf(d.pos)).attr("y", yOf(d.pos) - 16 * k)
        .attr("text-anchor", "middle").attr("font-size", Math.max(12 * k, 9))
        .attr("fill", C.gDark)
        .text(modText(d.mod));
    });
    letters.raise();

    /* cleavage polylines: prefix ions flag the row from above, suffix from below */
    const groups = d3.group(frags, f => f.ion + f.number);
    const marks = g.append("g").attr("fill", "none")
      .attr("stroke-linecap", "round").attr("stroke-linejoin", "round");
    for (const fs of groups.values()) {
      const frag = fs.reduce((a, b) => (b.intensity > a.intensity ? b : a));
      const idx = cleavageOf(frag, L);
      if (idx <= 0 || idx >= L) continue;
      if (rowOf(idx) !== rowOf(idx + 1)) continue;                  // row break
      if (byPos.get(idx)?.mod || byPos.get(idx + 1)?.mod) continue;  // mod adjacency
      const x = (xOf(idx) + xOf(idx + 1)) / 2, y = yOf(idx);
      const up = frag.ion === "a" || frag.ion === "b" || frag.ion === "c";
      const pts = up
        ? `${x - 2 * k},${y - 12 * k} ${x + 4 * k},${y - 10 * k} ${x + 4 * k},${y + 10 * k}`
        : `${x + 4 * k},${y - 10 * k} ${x + 4 * k},${y + 10 * k} ${x + 10 * k},${y + 13 * k}`;
      const mark = marks.append("polyline")
        .attr("class", "cleavage").datum(fkey(frag))
        .attr("points", pts).attr("stroke", ionColor(frag.ion))
        .attr("stroke-width", 1.2).attr("opacity", 0.95).attr("data-key", fkey(frag));
      mark.append("title").text(fkey(frag));
      /* generous invisible twin for hover precision */
      marks.append("polyline")
        .attr("class", "cleavage-hit")
        .attr("points", pts).attr("stroke", "transparent").attr("stroke-width", 9)
        .style("cursor", "pointer")
        .on("mouseenter", e => { o.onCleavage && o.onCleavage(e, frag); })
        .on("mouseleave", () => o.onCleavageLeave && o.onCleavageLeave());
    }

    return { g, height, xOf, yOf, marks: marks.selectAll("polyline.cleavage"), letters };
  }

  /* ---- ppm error panel: matched-fragment error vs observed m/z ----
   * o: {top, height, x (shared m/z scale), tol, labelEvery}
   * Domain symmetric, floor = stated tolerance so the tolerance frame is always
   * visible (the tolerance IS a claim — render it). */
  function errorPanel(svg, frags, o) {
    const C = forge.C;
    const errors = frags.map(f => f.ppm);
    const lim = Math.max(o.tol, Math.ceil(d3.max(errors.map(Math.abs)) / 10) * 10);
    const yOf = v => o.top + 16 + (1 - (v + lim) / (2 * lim)) * (o.height - 44);

    const g = svg.append("g");
    /* tolerance frame + zero line */
    [[-o.tol, C.threshold], [0, C.gFaint], [o.tol, C.threshold]].forEach(([v, col]) =>
      g.append("line")
        .attr("x1", o.x.range()[0]).attr("x2", o.x.range()[1])
        .attr("y1", yOf(v)).attr("y2", yOf(v))
        .attr("stroke", col).attr("stroke-width", v === 0 ? 1 : 0.8)
        .attr("stroke-dasharray", v === 0 ? "5 3" : "4 3"));
    [[-lim, "start"], [0, "middle"], [lim, "end"]].forEach(([v, anchor]) =>
      g.append("text").attr("x", o.x.range()[0] - 8).attr("y", yOf(v) + 3.5)
        .attr("text-anchor", anchor).attr("font-size", 10).attr("fill", C.gMid)
        .text(v > 0 ? "+" + v : v));
    g.append("text").attr("x", o.x.range()[0] - 8).attr("y", o.top + 12)
      .attr("text-anchor", "end").attr("font-size", 10).attr("fill", C.gSoft).text("ppm");

    const dots = g.selectAll(null).data(frags, fkey).join("circle")
      .attr("cx", f => o.x(f.mz)).attr("cy", f => yOf(f.ppm))
      .attr("r", 2.6).attr("fill", f => ionColor(f.ion)).attr("opacity", 0.85)
      .style("cursor", "pointer");
    dots.append("title").text(f => fragmentLabel(f));

    g.append("g").attr("transform", `translate(0,${o.top + o.height - 16})`)
      .call(d3.axisBottom(o.x).ticks(8).tickSize(3))
      .call(s => s.select(".domain").attr("stroke-width", C.swAxis));
    g.append("text")
      .attr("x", (o.x.range()[0] + o.x.range()[1]) / 2).attr("y", o.top + o.height + 2)
      .attr("text-anchor", "middle").attr("font-size", 11).attr("fill", C.gMid).text("m/z");

    return { g, dots, yOf };
  }

  /* ---- evidence annotations — the standard `annotations` contract field any
   * spectrum figure accepts. Two kinds, two layers; the caller owns layering:
   * bands go to a group appended BEFORE the stems, callouts AFTER the labels.
   *   {kind:"band", from, to, label?, color?}   tinted m/z region behind stems
   *   {kind:"callout", mz, label?, color?, dash?}  dashed flag at one m/z
   * Both redraw wholesale (the house clear-and-redraw rule); they take the
   * panel geometry so any spectrum template can mount them without changes. */
  function annBands(g, anns, o) {
    /* o: {x, top, height, fs} */
    const C = forge.C;
    g.selectAll("*").remove();
    anns.filter(a => a.kind === "band").forEach(a => {
      const x0 = o.x(a.from), x1 = o.x(a.to);
      g.append("rect")
        .attr("x", Math.min(x0, x1)).attr("y", o.top + 4)
        .attr("width", Math.max(Math.abs(x1 - x0), 1)).attr("height", o.height - 28)
        .attr("fill", a.color || C.gMid).attr("fill-opacity", 0.12);
      if (a.label) g.append("text")
        .attr("x", (x0 + x1) / 2).attr("y", o.top + 18)
        .attr("text-anchor", "middle").attr("font-size", o.fs || 10.5)
        .attr("fill", C.gDark).text(a.label);
    });
    return g;
  }

  function annCallouts(g, anns, o) {
    /* o: {x, top, baseline, fs} */
    const C = forge.C;
    g.selectAll("*").remove();
    anns.filter(a => a.kind === "callout").forEach(a => {
      const cx = o.x(a.mz);
      const line = g.append("line")
        .attr("x1", cx).attr("x2", cx)
        .attr("y1", o.top + 4).attr("y2", o.baseline)
        .attr("stroke", a.color || C.accent).attr("stroke-width", 1);
      if (a.dash !== false) line.attr("stroke-dasharray", "4,3");
      if (a.label) g.append("text")
        .attr("x", cx + 4).attr("y", o.top + 16)
        .attr("font-size", o.fs || 10.5).attr("font-weight", 700)
        .attr("fill", a.color || C.accent).text(a.label);
    });
    return g;
  }

  /* ---- spectrum panel bodies (promoted: two consumers — fragment-annotation,
   * spectrum-viewer). All take already-filtered items; colour defaults to the
   * ion tokens; identity via fkey so zoom re-joins stay cheap. ---- */
  function stems(g, items, o) {
    /* o: {x, y(f)->px tip, y0 px baseline, colorOf=ionColor, sw=1.6, key=fkey} */
    return g.selectAll("line.peak").data(items, o.key || fkey).join("line")
      .attr("class", "peak")
      .attr("x1", f => o.x(f.mz)).attr("x2", f => o.x(f.mz))
      .attr("y1", o.y0).attr("y2", f => o.y(f))
      .attr("stroke", f => (o.colorOf || ionColor)(f))
      .attr("stroke-width", o.sw ?? 1.6);
  }

  function peakHits(g, items, o) {
    /* o: {x, top, height, onMove(e,f), onOut, w=3, key=fkey} */
    return g.selectAll("rect.hit").data(items, o.key || fkey).join("rect")
      .attr("class", "hit")
      .attr("x", f => o.x(f.mz) - (o.w ?? 3)).attr("width", (o.w ?? 3) * 2)
      .attr("y", o.top).attr("height", o.height)
      .attr("fill", "transparent").style("cursor", "crosshair")
      .on("mousemove", (e, f) => o.onMove(e, f))
      .on("mouseout", () => o.onOut());
  }

  function topLabels(g, items, o) {
    /* o: {x, yTip(f)->px, rank=f.intensity, n=12, gap=34, dy=13, fs=10,
     *    text=f=>f.label, colorOf=ionColor} — top-n by rank, x-sorted, greedy
     *    2-level dodge (forge.dodgeX); white halo keeps labels legible on stems */
    const top = items.slice().sort((a, b) => (o.rank ? o.rank(b) - o.rank(a) : b.intensity - a.intensity))
      .slice(0, o.n ?? 12).sort((a, b) => o.x(a.mz) - o.x(b.mz));
    const levels = forge.dodgeX(top.map(f => o.x(f.mz)), o.gap ?? 34);
    const placed = [];
    top.forEach((f, i) => { if (levels[i] >= 0) placed.push({ f, lv: levels[i] }); });
    return g.selectAll(null).data(placed, d => fkey(d.f)).join("text")
      .attr("x", d => o.x(d.f.mz))
      .attr("y", d => o.yTip(d.f) - 8 - d.lv * (o.dy ?? 13))
      .attr("text-anchor", "middle").attr("font-size", o.fs ?? 10).attr("font-weight", 600)
      .attr("fill", d => (o.colorOf || ionColor)(d.f))
      .attr("paint-order", "stroke").attr("stroke", "#ffffff").attr("stroke-width", 2.4)
      .text(d => (o.text || (f => f.label))(d.f));
  }

  return { T, ionColor, fragmentLabel, fkey, cleavageOf, prsmTitle, errorPanel,
           annBands, annCallouts, stems, peakHits, topLabels };
})();

