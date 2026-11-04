// Envelope comparison — the ONE isotope-envelope figure for the runtime.
// Convention: observed peaks are grey sticks (never a connecting
// line); theoretical envelopes are colored lines with markers, one series per
// candidate, all on the same shared m/z axis. Apex registration lines (observed
// dashed dark, candidate dotted in its series color) make isotope-offset or
// calibration errors readable at a glance; the stats block (top-right) carries
// the per-candidate fit numbers. Data contract:
// { data: { panels: [{ title?, observed: {x[], y[]}, apex?: {observed, series[]},
//            series: [{label, color?, x[], y[], stat?}] }], xlabel, ylabel } }
// All extras are optional — a bare {observed, series} panel still renders.
// Colors come in the contract (assigned by the adapter); missing colors fall
// back to the tokens categorical palette. Hover shows the nearest point;
// drag = zoom, double-click = reset (house convention).
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft, zoomBrush } = forge;
  const D = window.DATA.data;
  const panels = D.panels || [];
  const xlabel = D.xlabel || "m/z";
  const ylabel = D.ylabel || "Relative intensity";

  const PANEL_H = 196, GAP = 70, TOP = 66, BOTTOM = 66, LEFT = 78;
  const F = forge.fig({ H: TOP + panels.length * PANEL_H + (panels.length - 1) * GAP + BOTTOM });
  const { W, H, k } = F;
  const R = W - 30 * k;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);
  const panelId = (svg.node().parentElement.id || "").replace(/^panel-/, "");

  const cat = i => C.cat[i % C.cat.length];
  const color = (s, i) => s.color || cat(i);

  // one shared x across panels (like the example figure)
  const allX = panels.flatMap(p => [p.observed.x, ...p.series.map(s => s.x)].flat());
  const xFull = d3.extent(allX);
  const x = d3.scaleLinear().domain(xFull).nice().range([LEFT, R]);
  const xHome = x.domain();   // nice-ified full view — the double-click reset target

  /* ---------- legend (top-left): observed + every series of panel 0 ---------- */
  const first = panels[0] || { series: [] };
  const entries = [{ label: "Observed (sticks)", color: C.gDark, dash: true },
    ...first.series.map(s => ({ label: s.label, color: s.color || cat(first.series.indexOf(s)) }))];
  let lx = LEFT;
  entries.forEach(e => {
    const g = svg.append("g");
    if (e.dash) g.append("line").attr("x1", lx).attr("x2", lx + 18 * k)
      .attr("y1", TOP - 38 * k).attr("y2", TOP - 38 * k)
      .attr("stroke", e.color).attr("stroke-width", 2);
    else g.append("circle").attr("cx", lx + 9 * k).attr("cy", TOP - 38 * k)
      .attr("r", 3.5 * k).attr("fill", e.color);
    g.append("text").attr("x", lx + 24 * k).attr("y", TOP - 34 * k)
      .attr("font-size", 10.5 * k).attr("fill", C.gDark).text(e.label);
    lx += 30 * k + e.label.length * 5.8 * k + 22 * k;
  });

  /* ---------- stats block (top-right, publication-style box) ---------- */
  const statLines = (first.series || [])
    .map((s, i) => ({ stat: s.stat, label: s.label, color: color(s, i) }))
    .filter(s => s.stat);
  if (statLines.length) {
    const box = svg.append("g");
    const fs = 10 * k, rowH = fs + 5 * k;
    const boxW = Math.max(...statLines.map(s =>
      (s.label.length + s.stat.length) * fs * 0.52)) + 20 * k;
    const boxH = statLines.length * rowH + 10 * k;
    const bx = R - boxW, by = TOP - boxH - 8 * k;
    box.append("rect").attr("x", bx).attr("y", by).attr("width", boxW)
      .attr("height", boxH).attr("fill", "#fff")
      .attr("stroke", C.gSoft).attr("rx", 3);
    statLines.forEach((s, i) => {
      const y = by + 8 * k + i * rowH + fs;
      box.append("circle").attr("cx", bx + 8 * k).attr("cy", y - fs * 0.32)
        .attr("r", 3 * k).attr("fill", s.color);
      box.append("text").attr("x", bx + 15 * k).attr("y", y)
        .attr("font-size", fs).attr("fill", C.ink)
        .text(`${s.label}  ${s.stat}`);
    });
  }

  /* ---------- panels: redrawable so brush-zoom + double-click reset work ---------- */
  function draw() {
    svg.select("g.panels").remove();
    const root = svg.append("g").attr("class", "panels");
    panels.forEach((p, pi) => {
      const yTop = TOP + pi * (PANEL_H + GAP);
      const yBot = yTop + PANEL_H;
      const allY = [p.observed.y, ...p.series.map(s => s.y)].flat().filter(Number.isFinite);
      const y = d3.scaleLinear().domain([0, (d3.max(allY) || 1) * 1.12]).nice()
        .range([yBot, yTop]);
      const g = root.append("g");

      // multi-panel ladders carry per-panel titles; single-panel pages already
      // have the page title from base.html — drawing it twice is noise
      if (panels.length > 1 && p.title)
        svg.append("text").attr("x", LEFT).attr("y", yTop - 14 * k)
          .attr("font-size", 12 * k).attr("font-weight", 700).attr("fill", C.ink)
          .text(p.title);

      /* apex registration lines: the absolute-position discriminator */
      if (p.apex && p.apex.observed != null)
        g.append("line").attr("x1", x(p.apex.observed)).attr("x2", x(p.apex.observed))
          .attr("y1", yTop).attr("y2", yBot)
          .attr("stroke", C.gDark).attr("stroke-width", 1.2).attr("stroke-dasharray", "5,3")
          .style("cursor", "help")
          .on("mousemove", e => showTip(tt, e, `<b>Observed apex</b><br>${xlabel}: ${p.apex.observed.toFixed(4)}`, 190))
          .on("mouseleave", () => hideTip(tt));
      (p.apex && p.apex.series || []).forEach((ax, si) => {
        if (ax == null) return;
        const s = p.series[si] || {};
        g.append("line").attr("x1", x(ax)).attr("x2", x(ax))
          .attr("y1", yTop).attr("y2", yBot)
          .attr("stroke", color(s, si)).attr("stroke-width", 1.1).attr("stroke-dasharray", "2,3")
          .style("cursor", "help")
          .on("mousemove", e => showTip(tt, e, `<b>${s.label || "candidate"} apex</b><br>${xlabel}: ${ax.toFixed(4)}`, 210))
          .on("mouseleave", () => hideTip(tt));
      });

      // observed: grey sticks (literature convention — never a connecting line)
      p.observed.x.forEach((mz, i) => {
        g.append("line").attr("x1", x(mz)).attr("x2", x(mz))
          .attr("y1", y(0)).attr("y2", y(p.observed.y[i]))
          .attr("stroke", C.gDark).attr("stroke-width", 1.6)
          .on("mousemove", e => showTip(tt, e,
            `<b>Observed</b><br>${xlabel}: ${mz}<br>${ylabel}: ${p.observed.y[i].toFixed(4)}`, 200))
          .on("mouseleave", () => hideTip(tt));
      });

      // theory series: line + haloed markers, one per candidate
      p.series.forEach((s, si) => {
        const pts = d3.zip(s.x, s.y);
        const ptsXY = pts.map(d => [x(d[0]), y(d[1])]);
        const selectFamily = () => {
          if (window.reportSelect && s.familyId) {
            window.reportSelect(panelId, {
              kind: "family", familyId: s.familyId,
              label: s.label, featureIds: s.featureIds || []
            });
          }
        };
        g.append("path").attr("d", d3.line()(ptsXY)).attr("fill", "none")
          .attr("stroke", color(s, si)).attr("stroke-width", 2)
          .style("cursor", s.familyId ? "pointer" : null)
          .on("click", selectFamily);
        g.selectAll(null).data(pts).join("circle")
          .attr("cx", d => x(d[0])).attr("cy", d => y(d[1])).attr("r", 3.2 * k)
          .attr("fill", color(s, si)).attr("stroke", "#fff").attr("stroke-width", 1)
          .style("cursor", s.familyId ? "pointer" : null)
          .on("click", selectFamily)
          .on("mousemove", (e, d) => showTip(tt, e,
            `<b>${s.label}</b>${s.stat ? `<br><span style="color:${C.gMid}">${s.stat}</span>` : ""}` +
            `<br>${xlabel}: ${d[0].toFixed(4)}<br>${ylabel}: ${d[1].toFixed(4)}`, 240))
          .on("mouseleave", () => hideTip(tt));
      });

      axLeft(root, y, { l: LEFT, n: 4 });
      axBottom(root, x, { y: yBot, n: 8, label: pi === panels.length - 1 ? xlabel : null });
    });
  }
  draw();

  svg.append("text").attr("transform", `rotate(-90)`)
    .attr("x", -(TOP + H / 2 - BOTTOM / 2)).attr("y", 18 * k)
    .attr("text-anchor", "middle").attr("font-size", 12 * k)
    .attr("fill", C.gMid).text(ylabel);

  const zb = zoomBrush(svg, {
    extent: [[LEFT, TOP], [R, H - BOTTOM]],
    x,
    onZoom: () => draw(),
  });
  svg.on("dblclick", () => zb.reset(xHome));   // house convention: double-click = reset
  window.__chartReady = true;
})();
