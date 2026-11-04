// Attention projection — horizontal bars, one panel per dataset, shared score
// scale (cross-dataset magnitudes stay honest). Bars are questions for review,
// not chemical identities. Top-ranked row takes the house accent; the rest sit
// in link blue with the score at the bar end. Contract:
// { data: { panels: [{ title, rows: [{ label, score, hover }] }], xlabel }, meta }
(function () {
  const { C, tip, showTip, hideTip, axBottom, gridV } = forge;
  const D = window.DATA.data;
  const panels = D.panels || [];
  const xlabel = D.xlabel || "attention score";

  const PANEL_W = 460, GAP = 70, TOP = 64, BOTTOM = 66, ROW_H = 26, BAR_H = 14;
  const F = forge.fig({ H: TOP + BOTTOM + Math.max(...panels.map(p => p.rows.length), 1) * ROW_H });
  const { W, H, k } = F;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip(svg);

  const xMax = Math.max(...panels.flatMap(p => p.rows.map(r => r.score)), 0.001) * 1.05;

  panels.forEach((p, pi) => {
    const colW = Math.min(PANEL_W * k, (W - 100 * k) / panels.length - GAP * k / panels.length);
    const cx = 96 * k + pi * (colW + GAP * k);
    const yTop = TOP + pi * 8 * k;
    const y = i => yTop + i * ROW_H * k;
    const x = d3.scaleLinear().domain([0, xMax]).range([cx + 4 * k, cx + colW]);
    const plotH = p.rows.length * ROW_H * k;

    svg.append("text").attr("x", cx).attr("y", yTop - 14 * k)
      .attr("font-size", 12 * k).attr("font-weight", 700).attr("fill", C.ink).text(p.title);
    gridV(svg.append("g"), x, { h: plotH, t: yTop, n: 5 });

    const g = svg.append("g").attr("class", "plotg");
    const top = Math.max(...p.rows.map(r => r.score), -Infinity);
    p.rows.forEach((r, i) => {
      const isTop = r.score >= top - 1e-9;
      const row = g.append("g").style("cursor", "crosshair")
        .on("mousemove", e => showTip(tt, e, `<b>${r.label}</b><br>${xlabel}: ${r.score.toFixed(3)}<br>${r.hover || ""}`, 250))
        .on("mouseleave", () => hideTip(tt));
      row.append("text").attr("x", cx - 8 * k).attr("y", y(i) + 15 * k)
        .attr("text-anchor", "end").attr("font-size", 11 * k)
        .attr("font-family", C.fontMono).attr("fill", C.gDark).text(r.label);
      row.append("rect").attr("x", x(0)).attr("y", y(i) + 4 * k)
        .attr("width", Math.max(x(r.score) - x(0), 1.5)).attr("height", BAR_H * k)
        .attr("rx", 2).attr("fill", isTop ? C.accent : C.link)
        .attr("fill-opacity", isTop ? 0.95 : 0.8);
      row.append("text").attr("x", x(r.score) + 7 * k).attr("y", y(i) + 15 * k)
        .attr("font-size", 10 * k).attr("font-family", C.fontMono)
        .attr("fill", C.gMid).text(r.score.toFixed(3));
    });

    axBottom(g, x, { y: yTop + plotH + 4 * k, n: 5,
                     label: pi === panels.length - 1 ? xlabel : null });
  });
  window.__chartReady = true;
})();
