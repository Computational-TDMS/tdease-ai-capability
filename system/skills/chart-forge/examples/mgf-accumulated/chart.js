// MGF accumulated-spectrum viewer — the quantify.py _sum_peak_lists contract
// running live: checked scans are chain-merged at tolerancePpm (intensity-
// weighted m/z, summed intensity), so ETD replicates sum into one boosted
// spectrum while the HCD series stays a toggle away. Grey = single scans
// (display-capped at the 400 strongest peaks; accumulation uses all peaks),
// accent = accumulated.
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, fmt, zoomBrush, legend } = forge;
  const D = window.DATA.data;
  const scans = D.scans;
  const PPM = D.tolerancePpm;
  const ANNS = D.annotations || [];
  const STICK_CAP = 400;

  const F = forge.fig({ H: 660 });
  const { W, H, k } = F;
  const SB = 712 * k;                       // sidebar x
  const M = { l: 64 * k, t: 30 * k };
  const PLOT_R = SB - 24 * k;               // spectrum panel right edge
  const yBot = H - 46 * k;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  /* ---------- accumulate: JS port of _sum_peak_lists (data-shaping, not chemistry).
     Chain-merge ascending peaks: within ppm of the last cluster → intensity-
     weighted m/z, summed intensity; otherwise a new cluster. ---------- */
  function accumulate(idxs) {
    const pts = [];
    idxs.forEach(i => scans[i].peaks.forEach(p => pts.push([p[0], p[1], i])));
    pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const cl = [];
    for (const [mz, inten, si] of pts) {
      const last = cl[cl.length - 1];
      if (!last || Math.abs(last.mz - mz) > mz * PPM / 1e6) {
        cl.push({ mz, int: inten, n: new Set([si]) });
      } else {
        const total = last.int + inten;
        last.mz = (last.mz * last.int + mz * inten) / total;
        last.int = total;
        last.n.add(si);
      }
    }
    return cl;
  }

  /* ---------- selection state: per-scan flags, default = the largest group */
  const groups = [];
  scans.forEach((s, i) => {
    let g = groups.find(g => g.name === s.activation);
    if (!g) groups.push(g = { name: s.activation, items: [] });
    g.items.push(i);
  });
  const nBiggest = Math.max(...groups.map(g => g.items.length));
  const scanOn = scans.map(s =>
    groups.find(g => g.name === s.activation).items.length === nBiggest);
  const selectedIdxs = () => scans.map((s, i) => i).filter(i => scanOn[i]);

  /* ---------- scales + axes */
  const full = d3.extent(scans.flatMap(s => s.peaks.map(p => p[0])));
  const x = d3.scaleLinear().domain(full).range([M.l, PLOT_R]);
  const y = d3.scaleLinear().range([yBot, M.t]);
  const gPlot = svg.append("g");
  const axX = axBottom(svg, x, { y: yBot, n: 10, label: "m/z" });
  const axY = axLeft(svg, y, { l: M.l, n: 5, label: "intensity (accumulated)" });

  let clustersNow = [];   // hover target: current accumulated peaks
  let mark = null;        // table-row highlight marker

  function stickPeaks(s) {
    if (s.peaks.length <= STICK_CAP) return s.peaks;
    return s.peaks.slice().sort((a, b) => b[1] - a[1]).slice(0, STICK_CAP)
      .sort((a, b) => a[0] - b[0]);
  }

  /* ---------- render: everything derives from scanOn + x.domain() */
  function render() {
    const sel = selectedIdxs();
    clustersNow = accumulate(sel);
    const inWin = clustersNow.filter(c => c.mz >= x.domain()[0] && c.mz <= x.domain()[1]);
    y.domain([0, (d3.max(inWin, c => c.int) || 1) * 1.08]).nice();

    gPlot.selectAll("*").remove();
    /* annotations: bands under the stems, callouts over them (forgeMs2 atoms) */
    forgeMs2.annBands(gPlot, ANNS, { x, top: M.t, height: yBot - M.t, fs: 10.5 * k });
    sel.forEach(i => {
      const g = gPlot.append("g");
      stickPeaks(scans[i]).forEach(p => g.append("line")
        .attr("x1", x(p[0])).attr("x2", x(p[0]))
        .attr("y1", yBot).attr("y2", y(p[1]))
        .attr("stroke", C.gMid).attr("stroke-width", 1).attr("stroke-opacity", 0.8));
    });
    gPlot.append("g").selectAll("line").data(clustersNow).join("line")
      .attr("x1", d => x(d.mz)).attr("x2", d => x(d.mz))
      .attr("y1", yBot).attr("y2", d => y(d.int))
      .attr("stroke", C.accent).attr("stroke-width", 1.6);
    forgeMs2.annCallouts(gPlot, ANNS, { x, top: M.t, baseline: yBot, fs: 10.5 * k });
    if (!sel.length) gPlot.append("text")
      .attr("x", (M.l + PLOT_R) / 2).attr("y", (M.t + yBot) / 2)
      .attr("text-anchor", "middle").attr("font-size", 13 * k).attr("fill", C.gMid)
      .text("no scans selected — tick a group on the right");

    axX.call(d3.axisBottom(x).ticks(10)).selectAll("text").attr("font-size", 11);
    axY.call(d3.axisLeft(y).ticks(4)).selectAll("text").attr("font-size", 11);
    drawPanel();
    drawTable(inWin, sel.length);
  }

  /* ---------- sidebar: precursor facts + activation-group checkboxes */
  function drawPanel() {
    const g = svg.select(".panel-side");
    g.selectAll("*").remove();
    const P = D.precursor;
    let ry = M.t + 10 * k;
    g.append("text").attr("x", SB).attr("y", ry).attr("font-size", 12.5 * k)
      .attr("font-weight", 700).attr("fill", C.ink)
      .text(`precursor ${P.mz.toFixed(4)}  ·  z = ${Math.abs(P.charge)}${P.charge >= 0 ? "+" : "-"}`);
    ry += 20 * k;
    g.append("text").attr("x", SB).attr("y", ry).attr("font-size", 11 * k).attr("fill", C.gMid)
      .text(`merge ${PPM} ppm chain · intensity-weighted m/z`);

    groups.forEach(gr => {
      ry += 26 * k;
      const allOn = gr.items.every(i => scanOn[i]);
      const gh = g.append("g").style("cursor", "pointer")
        .on("click", () => { gr.items.forEach(i => scanOn[i] = !allOn); render(); });
      gh.append("rect").attr("x", SB).attr("y", ry - 11 * k).attr("width", 10 * k).attr("height", 10 * k)
        .attr("fill", allOn ? C.accent : "none").attr("stroke", C.gMid).attr("stroke-width", 1);
      gh.append("text").attr("x", SB + 18 * k).attr("y", ry).attr("font-size", 12 * k)
        .attr("font-weight", 700).attr("fill", allOn ? C.ink : C.gMid)
        .text(`${gr.name} · ${gr.items.length}`);
      gr.items.forEach(i => {
        ry += 19 * k;
        const s = scans[i];
        const rh = g.append("g").style("cursor", "pointer")
          .on("click", () => { scanOn[i] = !scanOn[i]; render(); });
        rh.append("rect").attr("x", SB + 14 * k).attr("y", ry - 9 * k).attr("width", 8 * k).attr("height", 8 * k)
          .attr("fill", scanOn[i] ? C.gDark : "none").attr("stroke", C.gMid).attr("stroke-width", 0.8);
        rh.append("text").attr("x", SB + 30 * k).attr("y", ry).attr("font-size", 10.5 * k)
          .attr("font-family", C.fontMono).attr("fill", scanOn[i] ? C.gDark : C.gSoft)
          .text(`${s.scan} · ${s.rt.toFixed(1)} min · ${s.peaks.length} pk`);
      });
    });
    return ry;
  }

  /* ---------- top-peaks table (current x window, current selection) */
  function drawTable(inWin, nSel) {
    let ty = M.t + 30 * k + 26 * k;
    groups.forEach(gr => { ty += 26 * k + gr.items.length * 19 * k; });
    ty += 16 * k;
    const g = svg.select(".panel-table");
    g.selectAll("*").remove();
    const top = inWin.slice().sort((a, b) => b.int - a.int).slice(0, 14);
    g.append("text").attr("x", SB).attr("y", ty).attr("font-size", 12 * k)
      .attr("font-weight", 700).attr("fill", C.ink)
      .text(`top accumulated peaks (${nSel} scans)`);
    ["m/z", "intensity", "scans"].forEach((h, c) => g.append("text")
      .attr("x", SB + c * 100 * k).attr("y", ty + 20 * k).attr("font-size", 10 * k)
      .attr("fill", C.gMid).text(h));
    top.forEach((c, r) => {
      const row = g.append("g").style("cursor", "crosshair")
        .on("mouseenter", () => highlight(c))
        .on("mouseleave", () => highlight(null));
      [c.mz.toFixed(4), fmt("si")(c.int), `${c.n.size}/${nSel}`].forEach((v, ci) => row.append("text")
        .attr("x", SB + ci * 100 * k).attr("y", ty + 38 * k + r * 16 * k).attr("font-size", 10 * k)
        .attr("font-family", C.fontMono).attr("fill", C.gDark).text(v));
    });
  }

  function highlight(c) {
    if (mark) { mark.remove(); mark = null; }
    if (!c) return;
    mark = svg.append("g");
    mark.append("line").attr("x1", x(c.mz)).attr("x2", x(c.mz)).attr("y1", M.t).attr("y2", yBot)
      .attr("stroke", C.gGhost).attr("stroke-dasharray", "3 3");
    mark.append("circle").attr("cx", x(c.mz)).attr("cy", y(c.int)).attr("r", 4)
      .attr("fill", "none").attr("stroke", C.ink).attr("stroke-width", 1.4);
  }

  /* ---------- hover: nearest accumulated peak */
  svg.append("rect")
    .attr("x", M.l).attr("y", M.t).attr("width", PLOT_R - M.l).attr("height", yBot - M.t)
    .attr("fill", "transparent")
    .on("mousemove", (e) => {
      const [mx] = d3.pointer(e);
      let best = null, bd = Infinity;
      for (const c of clustersNow) {
        const d = Math.abs(x(c.mz) - mx);
        if (d < bd) { bd = d; best = c; }
      }
      if (!best || bd > 6) { hideTip(tt); return; }
      showTip(tt, e,
        `<b>m/z ${best.mz.toFixed(4)}</b><br>${fmt("si")(best.int)} · ${best.n.size} scan${best.n.size > 1 ? "s" : ""} contribute`,
        220);
    })
    .on("mouseleave", () => hideTip(tt));

  /* ---------- legend + brush zoom + reset */
  legend(svg, [
    { label: "single scans (capped)", color: C.gSoft },
    { label: "accumulated", color: C.accent },
  ], { x: M.l, y: M.t - 14 * k, pitch: 190 * k });

  const zb = zoomBrush(svg, {
    extent: [[M.l, M.t], [PLOT_R, yBot]],
    x,
    onZoom: () => render(),
  });
  svg.on("dblclick", () => zb.reset(full));   // house convention: double-click = reset

  svg.append("g").attr("class", "panel-side");
  svg.append("g").attr("class", "panel-table");
  render();
  window.__chartReady = true;
})();
