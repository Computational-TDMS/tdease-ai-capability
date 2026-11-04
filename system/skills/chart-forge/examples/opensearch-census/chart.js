// tdp-opensearch figures.py census, reproduced — visual contract matched to the
// plotly original: linear PrSM-count axis with 1.18 headroom, uniform sticks,
// evidence marks as dots ON their cluster's stick top with top-center labels
// (deep 13px coloured, routine 10px grey), edge labels anchor inward.
(function () {
  const { C, tip, showTip, hideTip, axBottom, axLeft } = forge;
  const { clusters, marks, xlim, floor } = window.DATA.data;
  const F = forge.fig({ W: 1280, H: 247 });
  const { W, H, k } = F;
  const M = { l: 64 * k, r: 42 * k, t: 34 * k, b: 44 * k };

  const VERDICT = { green: "#218739", orange: "#D55E00", grey: null };
  const markOf = m => {
    const deep = m.role === "deep";
    return VERDICT[m.color] || (deep ? "#D55E00" : C.gSoft);
  };

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);
  const panelId = (svg.node().parentElement.id || "").replace(/^panel-/, "");

  // page title comes from base.html (meta.title) — no hardcoded chart title

  const x = d3.scaleLinear().domain(xlim).range([M.l, W - M.r]);
  const inview = clusters.filter(c => c.inview);
  const yMax = d3.max(inview, c => c.count);
  const y = d3.scaleLinear().domain([0, yMax * 1.18]).range([H - M.b, M.t]);

  /* ---------- observation sticks: uniform, width 1; automatic UniMod
   * attributions color the stick (accent) and live in the tooltip — never
   * as scattered dots or on-axis text ---------- */
  function draw() {
    svg.select("g.sticks").remove();
    const g = svg.append("g").attr("class", "sticks");
    g.selectAll("line").data(inview, c => c.mass).join("line")
      .attr("x1", c => x(c.mass)).attr("x2", c => x(c.mass))
      .attr("y1", y(0)).attr("y2", c => y(c.count))
      .attr("stroke", c => c.annotation ? C.accent : C.gMid)
      .attr("stroke-width", 1)
      .style("cursor", "crosshair")
      .on("mousemove", (e, c) => showTip(tt, e,
        `<b>${c.mass > 0 ? "+" : ""}${c.mass.toFixed(6)} Da</b> · ${c.count} PrSM observations` +
        (c.samples ? ` · ${c.samples} samples` : "") +
        (c.localized != null ? ` · ${c.localized} localized` : "") +
        `<br>${c.dominant ? String(c.dominant).slice(0, 46) : ""}` +
        (c.residues && c.residues.length ? `<br>candidates: ${c.residues.join(", ")}` : "") +
        (c.best_e ? `<br>best E ${c.best_e.toExponential(1)}` : "") +
        (c.annotation ? `<br><b>${c.annotation}</b> <span style="color:${C.gSoft}">(automatic)</span>` : "") +
        (c.mark_labels.length ? `<br><b>${c.mark_labels.join(" · ")}</b> ${c.mark_verdicts.join("/")}` : ""),
        300))
      .on("click", (e, c) => {
        if (window.reportSelect) window.reportSelect(panelId, {
          kind: "mass-shift", mass: c.mass, count: c.count,
          assignmentIds: c.assignmentIds || []
        });
      })
      .on("mouseout", () => hideTip(tt));
    forge.legend(svg, [
      { label: "direct attribution", color: C.accent },
      { label: "no attribution", color: C.gMid },
    ], { x: M.l + 8, y: M.t - 14, pitch: 128, fs: 10 });
  }
  draw();

  /* ---------- evidence marks: dot on the stick top, label above.
   * Only curated deep/routine marks (explicit --marks input) ever draw;
   * automatic attributions are already shown by the stick color. ---------- */
  function drawMarks() {
    svg.select("g.marks").remove();
    const g = svg.append("g").attr("class", "marks");
    const span = xlim[1] - xlim[0];
    const occupied = [];
    const ordered = marks.filter(m => m.role === "deep" || m.role === "routine");
    ordered.forEach(m => {
      const cl = inview.find(c => Math.abs(c.mass - m.mass) < 0.01);
      if (!cl) return;
      const deep = m.role === "deep";
      const col = markOf(m);
      const px = x(m.mass), py = y(cl.count);
      g.append("circle")
        .attr("cx", px).attr("cy", py).attr("r", deep ? 5 : 4)
        .attr("fill", col)
        .attr("stroke", "#fff").attr("stroke-width", 1.2)
        .style("cursor", "help")
        .on("mousemove", e => showTip(tt, e,
          `<b>${m.label}</b> · ${m.role}${m.verdict ? " · " + m.verdict : ""}` +
          (m.observations != null ? `<br>${m.observations} cluster observations` : "") +
          `<br><span style="color:${C.gMid}">${m.basis}</span>`, 330))
        .on("mouseout", () => hideTip(tt));
      /* top-center label; near a frame edge the label anchors inward */
      const anchor = m.mass > xlim[1] - 0.12 * span ? "end"
        : m.mass < xlim[0] + 0.12 * span ? "start" : "middle";
      const lx = anchor === "start" ? px + 8 : anchor === "end" ? px - 8 : px;
      const fs = deep ? 13 : 10;
      const width = Math.max(22, m.label.length * fs * 0.56);
      const x0 = anchor === "start" ? lx : anchor === "end" ? lx - width : lx - width / 2;
      const x1 = x0 + width;
      let ly = py - (deep ? 12 : 9);
      while (occupied.some(box => x0 < box.x1 + 5 && x1 > box.x0 - 5 &&
             ly > box.y0 - 3 && ly - fs < box.y1 + 3)) {
        ly -= fs + 5;
      }
      occupied.push({x0, x1, y0: ly - fs, y1: ly});
      g.append("text")
        .attr("x", lx).attr("y", ly)
        .attr("text-anchor", anchor)
        .attr("font-size", fs)
        .attr("font-weight", deep ? 700 : 500)
        .attr("fill", deep ? col : C.gSoft)
        .attr("paint-order", "stroke").attr("stroke", "#fff").attr("stroke-width", 2.4)
        .text(m.label);
    });
  }
  drawMarks();

  axLeft(svg, y, { l: M.l, n: 6, label: "PrSM observations" });
  const axB = axBottom(svg, x, { y: H - M.b + 6 * k, n: 14, label: "Unexpected mass shift (Da)" });

  const zb = forge.zoomBrush(svg, {
    extent: [[M.l, M.t], [W - M.r, H - M.b]],
    x,
    onZoom: (dom) => {
      x.domain(dom);
      draw(); drawMarks();
      axB.call(d3.axisBottom(x).ticks(14));
      axB.select(".domain").attr("stroke-width", C.swAxis);
      axB.selectAll("text").attr("font-size", 11);
    },
  });
  svg.on("dblclick", () => zb.reset(xlim));   // house convention: double-click = reset
  window.__chartReady = true;
})();
