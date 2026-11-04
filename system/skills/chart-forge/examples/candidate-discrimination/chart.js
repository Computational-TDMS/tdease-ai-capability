/* One-flow candidate discrimination — single SVG, top to bottom:
 *   1. MS2 spectrum, stems coloured by which candidate claimed each peak
 *      (single = candidate colour, dark = claimed by both, grey = display noise);
 *   2. one block per candidate: its N-terminal series (c) error map above
 *      its OWN sequence strip, its C-terminal series (z) map below — every
 *      modification of the hypothesis is annotated on the strip (claims that
 *      differ between candidates in the candidate's colour), so each map's
 *      cliffs and gaps read against the marks that caused them;
 *   3. per-hypothesis ppm-error histograms — each candidate's own c/z (or
 *      b/y) overlaid in series colours, transparent, same direction.
 * One auditable object, N projections: peak <-> map cell <-> residue, across
 * every candidate at once — a missing counterpart is itself the evidence.
 * Data contract: ms2_visualize.py --mzml --candidates (kernel envelope matching). */
(function () {
  "use strict";
  const D = window.DATA.data;
  const C = forge.C;
  const CANDS = D.candidates, SEQ = D.sequence, TOL = D.tol, SCAN = D.scan;
  const N = SEQ.length, FOCUS = D.focusSite;
  const fsT = 11, fsL = 12.5;

  const ions = [...new Set(CANDS.flatMap(c => Object.keys(c.counts)))];
  const ionOrder = { c: 0, z: 1, b: 0, y: 1, a: 0, x: 1 };
  ions.sort((a, b) => (ionOrder[a] ?? 9) - (ionOrder[b] ?? 9) || a.localeCompare(b));

  const W = 1500;
  const M = { l: 76, r: 44 };

  const vikRamp = d3.interpolateRgbBasis(vik);
  const color = d3.scaleLinear().domain([-TOL, 0, TOL])
    .range([vikRamp(0), vikRamp(0.5), vikRamp(1)]).clamp(true);
  const tip = forge.tip();
  const ionColor = ion => ({ c: C.cIon, z: C.zIon, b: C.bIon, y: C.yIon,
    a: C.aIon, x: C.xIon })[ion] || C.accent;
  const SHARED_FILL = { ac: "#f6c344", me: "#f6c344", me2: "#f6c344", me3: "#f6c344",
    ox: "#8fd18f", cam: "#8fb8d1", su: "#d18fb8" };

  /* ---------- explicit layout ---------- */
  const specTop = 60, specH = 250, specBase = specTop + specH;   // 60..310
  const panelH = ion => (ion === "z" || ion === "y" || ion === "x" ? 100 : 92);
  const seqH = 54;
  // Per-candidate block: the N-terminal series map above the candidate's OWN
  // sequence strip, the C-terminal series map below it — the modifications
  // that define the hypothesis sit between its own maps, so the c-cliff and
  // z-gaps read directly against the strip's mod marks.
  const nTermIons = ions.filter(i => ["a", "b", "c"].includes(i));
  const cTermIons = ions.filter(i => !["a", "b", "c"].includes(i));
  const panelPos = {}, seqTopOf = {}, headerYOf = {};   // cid -> geometry
  let cursor = specBase + 34;
  CANDS.forEach(cand => {
    headerYOf[cand.id] = cursor + 14;
    let y = headerYOf[cand.id] + 12;
    panelPos[cand.id] = {};
    nTermIons.forEach(ion => {
      panelPos[cand.id][ion] = { top: y, h: panelH(ion) };
      y += panelH(ion) + 8;
    });
    seqTopOf[cand.id] = y + 10;          // headroom for mod tags above residues
    y = seqTopOf[cand.id] + seqH + 6;
    cTermIons.forEach(ion => {
      panelPos[cand.id][ion] = { top: y, h: panelH(ion) };
      y += panelH(ion) + 8;
    });
    cursor = y + 26;
  });
  const firstPanelTop = Math.min(...CANDS.flatMap(c =>
    Object.values(panelPos[c.id]).map(p => p.top)));
  const lastPanelBottom = Math.max(...CANDS.flatMap(c =>
    Object.values(panelPos[c.id]).map(p => p.top + p.h)));
  const axisY = lastPanelBottom + 14;                    // shared position axis

  /* ---------- bottom block: one histogram row per candidate ---------- */
  const histH = 100, histGap = 34, histTop = axisY + 54;
  const H = histTop + 14 + (CANDS.length - 1) * (histH + histGap) + histH + 52;
  const svg = d3.select("#chart").append("svg")
    .attr("viewBox", "0 0 " + W + " " + H).attr("width", W).attr("height", H)
    .style("font-family", C.font);

  /* ---------- indices ---------- */
  const mzMap = {}, cellMap = {}, peakMap = {};
  CANDS.forEach(c => {
    mzMap[c.id] = new Map();
    c.matches.forEach(m => {
      m.key = c.id + "|" + m.ion + m.n + "@" + m.q;
      m.cid = c.id;
      mzMap[c.id].set(Math.round(m.mz * 1e4).toString(), m);
    });
    cellMap[c.id] = new Map();
    peakMap[c.id] = new Map();
  });

  /* ---------- shared scales ---------- */
  const xSeq = d3.scaleLinear().domain([0, N]).range([M.l, W - M.r]);
  const pw = (W - M.l - M.r) / N;

  /* ---------- header ----------
     Candidate and ion identities are carried by the per-hypothesis block
     headers, the histogram row titles and the map row labels; no chip
     legends are stacked here. */
  svg.append("text").attr("x", M.l).attr("y", 24)
    .attr("font-size", 14).attr("font-weight", 700).attr("fill", C.gDark)
    .text((D.figureTitle || "Candidate discrimination") + " \u00b7 scan " + SCAN
      + (D.activation ? " \u00b7 " + D.activation : "")
      + " \u00b7 envelope-matched \u00b1" + TOL + " ppm");
  const cbW = 240, cbX = W - M.r - cbW, cbY = 40, cbH = 9;
  const grad = svg.append("defs").append("linearGradient").attr("id", "vcb");
  d3.range(31).forEach(i => grad.append("stop")
    .attr("offset", i * (100 / 30) + "%").attr("stop-color", color(-TOL + 2 * TOL * i / 30)));
  svg.append("rect").attr("x", cbX).attr("y", cbY).attr("width", cbW).attr("height", cbH)
    .attr("fill", "url(#vcb)").attr("rx", 2);
  [[0, "-" + TOL], [cbW / 2, "0"], [cbW, "+" + TOL]].forEach(pair =>
    svg.append("text").attr("x", cbX + pair[0]).attr("y", cbY - 4)
      .attr("text-anchor", pair[0] === 0 ? "start" : pair[0] === cbW ? "end" : "middle")
      .attr("font-size", fsT - 1).attr("fill", C.gMid).text(pair[1]));
  svg.append("text").attr("x", cbX + cbW + 6).attr("y", cbY + 8)
    .attr("font-size", fsT).attr("fill", C.gMid).text("ppm");

  /* ---------- highlight machinery ---------- */
  const clLines = [];
  function clearHi() {
    svg.selectAll(".pk-m")
      .attr("stroke", function () { return d3.select(this).attr("data-c"); })
      .attr("stroke-width", function () { return d3.select(this).attr("data-w"); });
    svg.selectAll(".cell-m").attr("stroke", null).attr("stroke-width", null);
    svg.selectAll(".pk-un").attr("opacity", 1);
    svg.selectAll(".pk-m, .cell-m").attr("opacity", 1);
    clLines.forEach(l => l.remove()); clLines.length = 0;
  }
  function cleave(cl) {
    const x = xSeq(cl) + pw / 2;
    clLines.push(svg.append("line").attr("x1", x).attr("x2", x)
      .attr("y1", firstPanelTop - 6).attr("y2", axisY)
      .attr("stroke", C.gDark).attr("stroke-width", 1).attr("stroke-dasharray", "3,3")
      .attr("pointer-events", "none"));
  }
  function highlight(matches) {
    svg.selectAll(".pk-un").attr("opacity", 0.10);
    matches.forEach(m => {
      const pe = peakMap[m.cid].get(m.key), ce = cellMap[m.cid].get(m.key);
      if (pe) pe.attr("stroke", C.gDark).attr("stroke-width", 1.3);
      if (ce) ce.attr("stroke", C.gDark).attr("stroke-width", 1.4);
    });
    [...new Set(matches.map(m => m.cleav))].forEach(cleave);
  }
  function coverTest(m, r) {
    return ["c", "a", "b"].includes(m.ion) ? m.n >= r : m.cleav < r;
  }
  function hiResidue(r, cid) {
    svg.selectAll(".pk-un").attr("opacity", 0.06);
    CANDS.forEach(c => c.matches.forEach(m => {
      if (c.id !== cid) return;          // this strip speaks for ONE hypothesis
      const on = coverTest(m, r);
      const pe = peakMap[m.cid].get(m.key), ce = cellMap[m.cid].get(m.key);
      if (pe) pe.attr("opacity", on ? 1 : 0.12);
      if (ce) ce.attr("opacity", on ? 1 : 0.15);
    }));
    cleave(r - 0.2);
  }

  /* ---------- spectrum panel ---------- */
  const maxI = d3.max(D.peaks, p => p.i) || 1;
  const xMz = d3.scaleLinear()
    .domain([Math.floor(d3.min(D.peaks, p => p.mz)), Math.ceil(d3.max(D.peaks, p => p.mz))])
    .range([M.l, W - M.r]);
  const ySpec = d3.scaleLinear().domain([0, maxI]).range([specBase, specTop]);
  svg.append("text").attr("x", M.l).attr("y", specTop - 10)
    .attr("font-size", fsL).attr("font-weight", 700).attr("fill", C.gDark)
    .text("MS2 spectrum \u2014 grey unmatched \u00b7 coloured by ion series "
      + "(same encoding as the histograms below) \u00b7 dark = claimed by both");
  const claimOf = new Map();
  CANDS.forEach(c => c.matches.forEach(m => {
    const k = Math.round(m.mz * 1e4).toString();
    if (!claimOf.has(k)) claimOf.set(k, { set: new Set(), any: [] });
    const rec = claimOf.get(k);
    rec.set.add(c.id); rec.any.push(m);
  }));
  const gUn = svg.append("g");
  gUn.selectAll("line").data(D.peaks.filter(p => !claimOf.has(Math.round(p.mz * 1e4).toString())))
    .join("line").attr("class", "pk pk-un")
    .attr("x1", p => xMz(p.mz)).attr("x2", p => xMz(p.mz))
    .attr("y1", specBase).attr("y2", p => ySpec(p.i))
    .attr("stroke", "#c9c9c9").attr("stroke-width", 0.7);
  const gM = svg.append("g");
  claimOf.forEach((rec, k) => {
    const top = rec.any.reduce((a, m) => (m.i > a.i ? m : a), rec.any[0]);
    const line = gM.append("line").attr("class", "pk pk-m")
      .attr("data-c", ionColor(top.ion)).attr("data-w", rec.set.size > 1 ? 1.9 : 1.5)
      .attr("x1", xMz(top.mz)).attr("x2", xMz(top.mz))
      .attr("y1", specBase).attr("y2", ySpec(top.i))
      .attr("stroke", ionColor(top.ion)).attr("stroke-width", rec.set.size > 1 ? 1.9 : 1.5);
    rec.any.forEach(m => peakMap[m.cid].set(m.key, line));
    line.on("mouseenter", function (e) {
        highlight(rec.any);
        const lines = CANDS.map(c => {
          const m = rec.any.find(x => x.cid === c.id);
          return m ? "<span style='color:" + c.color + "'>\u25a0 " + c.label + "</span> \u2014 <b>"
            + m.ion + m.n + "</b> \u00b7 " + m.q + "+ \u00b7 "
            + (m.ppm >= 0 ? "+" : "") + m.ppm.toFixed(2) + " ppm"
            : "<span style='color:" + c.color + "'>\u25a1 " + c.label + "</span> \u2014 not matched";
        });
        forge.showTip(tip, e, "<b>m/z " + top.mz.toFixed(3) + "</b> \u00b7 i " + top.i.toFixed(0)
          + "<br>" + lines.join("<br>"));
      })
      .on("mouseleave", function () { clearHi(); forge.hideTip(tip); });
  });
  const seenMz = new Set();
  CANDS.flatMap(c => c.matches).sort((a, b) => b.i - a.i)
    .filter(m => {
      const k = Math.round(m.mz * 1e4).toString();
      if (seenMz.has(k)) return false;
      seenMz.add(k);
      return true;
    })
    .slice(0, 18)
    .forEach(m => {
      svg.append("text").attr("x", xMz(m.mz)).attr("y", ySpec(m.i) - 5)
        .attr("transform", "rotate(-62 " + xMz(m.mz) + " " + (ySpec(m.i) - 5) + ")")
        .attr("font-size", fsT - 1).attr("fill", ionColor(m.ion)).attr("font-weight", 700)
        .text(m.ion.toUpperCase() + m.n);
    });
  svg.append("g").attr("transform", "translate(0," + specBase + ")")
    .call(d3.axisBottom(xMz).ticks(9).tickSize(4))
    .call(g => g.select(".domain").attr("stroke", C.gFaint))
    .call(g => g.selectAll("text").attr("font-size", fsT).attr("fill", C.gDark));
  svg.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", specBase + 32)
    .attr("text-anchor", "middle").attr("font-size", fsL).attr("fill", C.gDark).text("m/z");
  /* Relative intensity, as a percentage of the base peak. Without an ordinate
     the panel cannot be read for abundance: a stem at 4 % and a stem at 40 %
     look alike in a short panel, and the reader has no way to tell a sparse
     region from a drawing fault. */
  svg.append("g").attr("transform", "translate(" + M.l + ",0)")
    .call(d3.axisLeft(ySpec)
      .tickValues([0, 0.25, 0.5, 0.75, 1].map(f => f * maxI))
      .tickSize(4).tickFormat(v => Math.round(v / maxI * 100)))
    .call(g => g.select(".domain").attr("stroke", C.gFaint))
    .call(g => g.selectAll("text").attr("font-size", fsT).attr("fill", C.gDark));
  svg.append("text").attr("transform", "rotate(-90)")
    .attr("x", -(specTop + specBase) / 2).attr("y", M.l - 58)
    .attr("text-anchor", "middle").attr("font-size", fsL).attr("fill", C.gDark)
    .text("relative intensity (% base peak)");

  /* ---------- error-map panels ---------- */
  function candHeader(c, y, headerTop) {
    const gaps = ions.map(ion => {
      const g = (c.worstGap || {})[ion];
      return g && g.length >= 3
        ? ion + " gap " + g.start + "\u2013" + g.end + " (" + g.length + ")"
        : null;
    }).filter(Boolean);
    const site = (c.siteIons || []).map(s =>
      s.ion + s.number + (s.states.length ? "\u2713" : "\u2717")).join("/");
    svg.append("text").attr("x", M.l).attr("y", y)
      .attr("font-size", fsL).attr("font-weight", 700).attr("fill", c.color)
      .text(c.label + " \u2014 " + c.total + " matched charge states (positions: "
        + ions.map(i => i + " " + (c.counts[i] || 0)).join(" \u00b7 ") + ")"
        + (gaps.length ? " \u00b7 " + gaps.join(" \u00b7 ") : " \u00b7 no \u22653 gap")
        + (site ? " \u00b7 site " + site : ""));
    return headerTop;
  }
  function mapPanel(c, ion) {
    const pos = panelPos[c.id][ion];
    const top = pos.top, h = pos.h;
    svg.append("text").attr("x", M.l - 8).attr("y", top - 4)
      .attr("text-anchor", "end").attr("font-size", fsT - 1)
      .attr("font-weight", 700).attr("fill", C.gMid).text(ion + " \u2192");
    const rows = D.charges;
    const ch = h / rows.length;
    const y = d3.scaleBand().domain(rows).range([top, top + h]);
    rows.forEach(q => svg.append("text").attr("x", M.l - 8).attr("y", y(q) + ch / 2 + 3.5)
      .attr("text-anchor", "end").attr("font-size", fsT - 1).attr("fill", C.gMid).text(q + "+"));
    rows.forEach(q => { for (let p = 0; p < N; p++) {
      svg.append("rect").attr("x", xSeq(p)).attr("y", y(q))
        .attr("width", Math.max(pw - 0.6, 1)).attr("height", Math.max(ch - 0.6, 1))
        .attr("fill", C.gap);
    }});
    c.matches.filter(m => m.ion === ion).forEach(m => {
      const rect = svg.append("rect").attr("class", "cell-m")
        .attr("x", xSeq(m.cleav)).attr("y", y(m.q))
        .attr("width", Math.max(pw - 0.6, 1)).attr("height", Math.max(ch - 0.6, 1))
        .attr("fill", color(m.ppm));
      cellMap[c.id].set(m.key, rect);
      rect.on("mouseenter", function (e) {
          highlight([m]);
          const others = CANDS.filter(x => x.id !== c.id).map(x => {
            const twin = x.matches.find(t => t.ion === m.ion && t.n === m.n);
            if (twin) return "<span style='color:" + x.color + "'>\u25a0 " + x.label
              + "</span> " + twin.ion + twin.n + " \u00b7 " + twin.q + "+ \u00b7 "
              + (twin.ppm >= 0 ? "+" : "") + twin.ppm.toFixed(2) + " ppm";
            const gap = (x.missingRuns[m.ion] || [])
              .find(r => m.n >= r.start && m.n <= r.end);
            return "<span style='color:" + x.color + "'>\u25a1 " + x.label + "</span> "
              + m.ion + m.n + " MISSING"
              + (gap ? " (inside gap " + gap.start + "\u2013" + gap.end + ")" : "");
          });
          forge.showTip(tip, e,
            "<span style='color:" + c.color + "'>\u25a0 " + c.label + "</span> <b>"
            + m.ion + m.n + "</b> \u00b7 " + m.q + "+ \u00b7 cleave " + m.cleav
            + "<br>error " + (m.ppm >= 0 ? "+" : "") + m.ppm.toFixed(2) + " ppm"
            + (others.length ? "<br>" + others.join("<br>") : ""));
        })
        .on("mouseleave", function () { clearHi(); forge.hideTip(tip); });
    });
  }
  CANDS.forEach(cand => {
    candHeader(cand, headerYOf[cand.id]);
    ions.forEach(ion => mapPanel(cand, ion));
  });

  /* ---------- per-candidate sequence strips (own mods between own maps) ---------- */
  CANDS.forEach(cand => {
    const stripTop = seqTopOf[cand.id];
    const mods = new Map((cand.mods || [])
      .map(m => [m.nTerm ? "nTerm" : m.cTerm ? "cTerm" : m.pos, m]));
    svg.append("text").attr("x", M.l - 8).attr("y", stripTop + 32)
      .attr("text-anchor", "end").attr("font-size", fsT).attr("font-weight", 700)
      .attr("fill", C.gMid).text("SEQ");
    const nTerm = mods.get("nTerm");
    if (nTerm) svg.append("text").attr("x", M.l - 8).attr("y", stripTop + 12)
      .attr("text-anchor", "end").attr("font-size", fsT)
      .attr("fill", nTerm.diff ? cand.color : "#8a6d1a")
      .attr("font-weight", nTerm.diff ? 700 : 400)
      .text("[" + nTerm.tag + "]-");
    const cTerm = mods.get("cTerm");
    if (cTerm) svg.append("text").attr("x", xSeq(N) + 8).attr("y", stripTop + 12)
      .attr("font-size", fsT).attr("fill", cTerm.diff ? cand.color : "#8a6d1a")
      .attr("font-weight", cTerm.diff ? 700 : 400).text("-[" + cTerm.tag + "]");
    for (let r = 1; r <= N; r++) {
      const x = xSeq(r - 1);
      const mod = mods.get(r);
      const g = svg.append("g").attr("class", "res-g").style("cursor", "crosshair");
      if (mod) {
        g.append("rect").attr("x", x + 0.6).attr("y", stripTop + 6)
          .attr("width", Math.max(pw - 1.2, 1)).attr("height", seqH - 12)
          .attr("fill", mod.diff ? cand.color : (SHARED_FILL[mod.tag] || "#b9a7e8"))
          .attr("opacity", 0.55).attr("rx", 2);
        g.append("text").attr("x", x + pw / 2).attr("y", stripTop - 2)
          .attr("text-anchor", "middle").attr("font-size", 9.5)
          .attr("fill", mod.diff ? cand.color : "#8a6d1a")
          .attr("font-weight", mod.diff ? 700 : 400).text(mod.tag);
      }
      const letter = g.append("text").attr("x", x + pw / 2).attr("y", stripTop + 32)
        .attr("text-anchor", "middle").attr("font-size", 12)
        .attr("font-weight", mod ? 700 : 400).attr("fill", C.gDark).text(SEQ[r - 1]);
      g.on("mouseenter", function (e) {
          letter.attr("font-weight", 700);
          hiResidue(r, cand.id);
          const per = CANDS.map(c => {
            const cnt = ions.map(ion =>
              ion + c.matches.filter(m => m.ion === ion && coverTest(m, r)).length).join(" ");
            return "<span style='color:" + c.color + "'>" + c.label + "</span> " + cnt;
          });
          const modHere = mod
            ? " \u00b7 " + (mod.diff
                ? "<b style='color:" + cand.color + "'>" + mod.tag + "</b> (claim)"
                : mod.tag)
            : "";
          const focusNote = FOCUS === r ? " \u2190 focus site" : "";
          forge.showTip(tip, e, "<b>residue " + r + " = " + SEQ[r - 1] + "</b>" + modHere
            + focusNote + "<br>" + per.join(" \u00b7 "));
        })
        .on("mouseleave", function () {
          letter.attr("font-weight", mod ? 700 : 400); clearHi(); forge.hideTip(tip);
        });
    }
  });

  /* ---------- shared position axis ---------- */
  svg.append("g").attr("transform", "translate(0," + axisY + ")")
    .call(d3.axisBottom(xSeq).ticks(10).tickSize(4))
    .call(g => g.select(".domain").attr("stroke", C.gFaint))
    .call(g => g.selectAll("text").attr("font-size", fsT).attr("fill", C.gDark));
  svg.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", axisY + 20)
    .attr("text-anchor", "middle").attr("font-size", fsT).attr("fill", C.gMid)
    .text("fragment position (cleavage site) \u00b7 c/a/b ions cleave after, x/y/z ions cover to the C-terminus");

  /* ---------- focus guide (behind everything) ----------
     The guide marks a residue, so it runs only where the position axis is the
     abscissa — the same span `cleave()` uses. Starting it at specTop put a
     sequence coordinate on the spectrum's m/z axis, where it read as an
     unlabelled m/z marker that the focus site has nothing to do with. */
  if (FOCUS) {
    svg.insert("line", ":first-child")
      .attr("x1", xSeq(FOCUS - 1) + pw / 2).attr("x2", xSeq(FOCUS - 1) + pw / 2)
      .attr("y1", firstPanelTop - 6).attr("y2", axisY)
      .attr("stroke", C.gDark).attr("stroke-width", 0.8)
      .attr("stroke-dasharray", "2,4").attr("opacity", 0.45)
      .attr("pointer-events", "none");
  }

  /* ---------- per-candidate ppm-error histograms (c/z overlaid) ---------- */
  // One row per hypothesis; its own c/z (or b/y) series are overlaid there —
  // transparent fills in the same series colours the spectrum uses, bars
  // growing the same direction from one baseline. The primary read is
  // within-row: a c/z mean divergence flags systematic disagreement with that
  // model (same signal as the kernel's seriesDivergence); rows stacked on one
  // ppm scale then compare candidates.
  svg.append("text").attr("x", M.l).attr("y", histTop - 10)
    .attr("font-size", fsL).attr("font-weight", 700).attr("fill", C.gDark)
    .text("per-hypothesis ppm-error histograms \u2014 its own c/z overlaid in "
      + "series colours \u00b7 a displaced series mean indicates systematic "
      + "disagreement with that model");
  const lastBase = histTop + 14 + (CANDS.length - 1) * (histH + histGap) + histH;
  CANDS.forEach((cand, ci) => {
    const top = histTop + 14 + ci * (histH + histGap), base = top + histH;
    const lo = -TOL * 1.15, hi = TOL * 1.15;
    const x = d3.scaleLinear().domain([lo, hi]).range([M.l, W - M.r]);
    const makeBins = d3.bin().value(m => m.ppm).domain([lo, hi]).thresholds(23);
    const perIon = ions.map(ion => {
      const series = cand.matches.filter(m => m.ion === ion);
      return { ion, series, bins: makeBins(series) };
    }).filter(p => p.series.length);
    const y = d3.scaleLinear()
      .domain([0, d3.max(perIon.flatMap(p => p.bins), b => b.length) || 4])
      .range([base, top]);

    const title = svg.append("text").attr("x", M.l).attr("y", top - 5)
      .attr("font-size", fsL - 1).attr("font-weight", 700);
    title.append("tspan").attr("fill", cand.color).text(cand.label);
    perIon.forEach(p => {
      const mean = d3.mean(p.series, m => m.ppm);
      title.append("tspan").attr("fill", C.gMid).text("   \u00b7  ");
      title.append("tspan").attr("fill", ionColor(p.ion))
        .text(p.ion + " \u03bc" + (mean >= 0 ? "+" : "") + mean.toFixed(1)
          + " n=" + p.series.length);
    });

    svg.append("line").attr("x1", x(0)).attr("x2", x(0))
      .attr("y1", top).attr("y2", base)
      .attr("stroke", C.gFaint).attr("stroke-width", 1);
    [-TOL, TOL].forEach(v => svg.append("line")
      .attr("x1", x(v)).attr("x2", x(v))
      .attr("y1", top).attr("y2", base)
      .attr("stroke", C.gFaint).attr("stroke-dasharray", "2,3").attr("opacity", 0.8));

    perIon.forEach(p => {
      if (p.series.length < 3) return;
      const mean = d3.mean(p.series, m => m.ppm);
      svg.append("line").attr("x1", x(mean)).attr("x2", x(mean))
        .attr("y1", top).attr("y2", base)
        .attr("stroke", ionColor(p.ion)).attr("stroke-width", 1.2)
        .attr("stroke-dasharray", "5,3").attr("opacity", 0.9)
        .style("pointer-events", "none");
    });

    perIon.forEach(p => p.bins.forEach(b => {
      if (!b.length) return;
      svg.append("rect")
        .attr("x", x(b.x0) + 0.75)
        .attr("y", y(b.length))
        .attr("width", Math.max(x(b.x1) - x(b.x0) - 1.5, 0.8))
        .attr("height", Math.max(base - y(b.length), 0.8))
        .attr("fill", ionColor(p.ion)).attr("opacity", 0.55)
        .attr("stroke", ionColor(p.ion)).attr("stroke-width", 0.8)
        .on("mouseenter", e => {
          highlight(b);
          const sameBin = perIon.filter(o => o !== p).map(o => {
            const n = o.bins.find(bb => bb.x0 === b.x0)?.length || 0;
            return "<span style='color:" + ionColor(o.ion) + "'>"
              + (n ? "\u25a0 " : "\u25a1 ") + o.ion + "-series</span> n=" + n;
          });
          forge.showTip(tip, e,
            "<span style='color:" + cand.color + "'>\u25a0 " + cand.label
            + "</span> <b>" + p.ion + "-series</b><br>error "
            + b.x0.toFixed(1) + " .. " + b.x1.toFixed(1) + " ppm \u00b7 n=" + b.length
            + (sameBin.length ? "<br>same bin \u00b7 " + sameBin.join(" \u00b7 ") : ""));
        })
        .on("mouseleave", () => { clearHi(); forge.hideTip(tip); });
    }));

    svg.append("g").attr("transform", "translate(0," + base + ")")
      .call(d3.axisBottom(x).ticks(5).tickSize(4))
      .call(g => g.select(".domain").attr("stroke", C.gFaint))
      .call(g => g.selectAll("text").attr("font-size", fsT).attr("fill", C.gDark));
    svg.append("g").attr("transform", "translate(" + M.l + ",0)")
      .call(d3.axisLeft(y).ticks(3).tickSize(3))
      .call(g => g.select(".domain").attr("stroke", C.gFaint))
      .call(g => g.selectAll("text").attr("font-size", fsT - 1).attr("fill", C.gMid));
  });
  svg.append("text").attr("x", (M.l + W - M.r) / 2).attr("y", lastBase + 32)
    .attr("text-anchor", "middle").attr("font-size", fsL).attr("fill", C.gDark)
    .text("ppm error");

  svg.append("text").attr("x", M.l).attr("y", H - 8)
    .attr("font-size", fsT).attr("fill", C.gSoft)
    .text(D.peaks.length + " displayed peaks \u00b7 "
      + (FOCUS ? "focus site @" + FOCUS : ""));

  window.__chartReady = true;
})();
