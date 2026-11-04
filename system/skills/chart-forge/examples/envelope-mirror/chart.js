// Envelope mirror — the default view for `envelope-regions` extractions.
// One card per source JSON (a candidate model, a scan, a charge state);
// several cards in one figure is the model-vs-model / scan-vs-scan comparison.
// Card anatomy: colored title + verdict subtitle; observed peaks up (black
// sticks, in-slot peaks carry their ±ppm label); theoretical comb down (the
// card color, M+k slot labels); symmetric ±100% axis so a displaced or
// shape-wrong envelope reads as tilt and gaps, never as a scale change.
// Data contract:
// { data: { panels: [{ title, subtitle?, color?,
//            observed: {x[], y[], ppm?[], matched?[]},
//            theoretical: {x[], y[], offsets[]} }], xlabel?, ylabel? } }
// y values are percent of the panel's own maximum on both halves; colors come
// in the contract (assigned by the adapter), falling back to the categorical
// tokens. Hover shows m/z, intensity and ppm; a ppm label marks an observed
// peak inside its theoretical slot (|err| <= adapter tolerance).
(function () {
  const { C, tip, showTip, hideTip, dodgeX } = forge;
  const D = window.DATA.data;
  const panels = D.panels || [];
  const xlabel = D.xlabel || "m/z";
  const ylabel = D.ylabel || "rel. intensity (%)";

  const F = forge.fig({ H: 10 });            // width + k from journal/tokens
  const { W, k } = F;
  const COLS = panels.length > 1 ? 2 : 1;
  const GAP = 18 * k, PAD_T = 10 * k;
  const CARD_W = (W - GAP * (COLS - 1)) / COLS;
  const CARD_H = 320 * k;
  const M_L = 60 * k, M_R = 14 * k, M_T = 56 * k, AXIS_H = 40 * k;
  const PLOT_H = CARD_H - M_T - AXIS_H;
  const rows = Math.ceil(panels.length / COLS);
  const H = PAD_T + rows * CARD_H + (rows - 1) * GAP + 12 * k;

  const svg = d3.select("#chart").append("svg")
    .attr("width", W).attr("height", H).attr("viewBox", `0 0 ${W} ${H}`);
  const tt = tip();
  const color = (p, i) => p.color || C.cat[i % C.cat.length];

  panels.forEach((p, i) => {
    const col = i % COLS, row = Math.floor(i / COLS);
    const cx = col * (CARD_W + GAP), cy = PAD_T + row * (CARD_H + GAP);
    const pc = color(p, i);
    const l = cx + M_L, r = cx + CARD_W - M_R;
    const top = cy + M_T, mid = top + PLOT_H / 2, bot = top + PLOT_H;

    const card = svg.append("g");
    card.append("rect").attr("x", cx).attr("y", cy)
      .attr("width", CARD_W).attr("height", CARD_H)
      .attr("rx", 10 * k).attr("fill", "#fff")
      .attr("stroke", C.gSoft).attr("stroke-width", 1);
    card.append("text").attr("x", cx + 16 * k).attr("y", cy + 26 * k)
      .attr("font-size", 13 * k).attr("font-weight", 700).attr("fill", pc)
      .text(p.title);
    if (p.subtitle)
      card.append("text").attr("x", cx + 16 * k).attr("y", cy + 43 * k)
        .attr("font-size", 11 * k).attr("fill", C.gMid).text(p.subtitle);

    // per-card x: union of both sides, padded so edge sticks keep their labels
    const xs = [...p.observed.x, ...p.theoretical.x];
    const [x0, x1] = d3.extent(xs);
    const padX = Math.max((x1 - x0) * 0.03, 1e-6);
    const x = d3.scaleLinear().domain([x0 - padX, x1 + padX]).range([l, r]);
    // both halves normalize to their own max = 100 but DRAW to a 125 domain:
    // the tallest stick on either side stops short of the card edge, so its
    // label (ppm above, M+k below) keeps clear space and never rides the axis
    const up = v => mid - v * (PLOT_H / 2) / 125;
    const down = v => mid + v * (PLOT_H / 2) / 125;
    const y = d3.scaleLinear().domain([-125, 125]).range([bot, top]);

    card.append("line").attr("x1", l).attr("x2", r)
      .attr("y1", mid).attr("y2", mid)
      .attr("stroke", C.gDark).attr("stroke-width", 1.2);
    card.append("text").attr("x", l + 8 * k).attr("y", top + 12 * k)
      .attr("font-size", 10 * k).attr("font-style", "italic").attr("fill", C.gMid)
      .text("observed");
    card.append("text").attr("x", l + 8 * k).attr("y", bot - 6 * k)
      .attr("font-size", 10 * k).attr("font-style", "italic").attr("fill", C.gMid)
      .text("theoretical");

    // observed up: black when inside a theoretical slot, gray otherwise
    const ppmOf = j => (p.observed.ppm || [])[j];
    const inSlot = j => p.observed.matched
      ? !!p.observed.matched[j]
      : Number.isFinite(ppmOf(j));
    p.observed.x.forEach((mz, j) => {
      const m = inSlot(j);
      card.append("line").attr("x1", x(mz)).attr("x2", x(mz))
        .attr("y1", mid).attr("y2", up(p.observed.y[j]))
        .attr("stroke", m ? C.gDark : C.gSoft)
        .attr("stroke-opacity", m ? 1 : 0.55).attr("stroke-width", 1.7)
        .on("mousemove", e => showTip(tt, e,
          `<b>Observed</b><br>${xlabel}: ${mz.toFixed(4)}<br>` +
          `${ylabel.split(" ")[0]}.: ${p.observed.y[j].toFixed(1)}%` +
          (Number.isFinite(ppmOf(j)) ? `<br>slot err ${ppmOf(j) >= 0 ? "+" : ""}${ppmOf(j).toFixed(1)} ppm` : ""), 210))
        .on("mouseleave", () => hideTip(tt));
    });
    // ±ppm labels above in-slot observed peaks (two dodge levels)
    const lab = p.observed.x.map((mz, j) => ({ mz, j, ppm: ppmOf(j), y: p.observed.y[j] }))
      .filter(o => Number.isFinite(o.ppm) && inSlot(o.j));
    const levels = dodgeX(lab.map(o => x(o.mz)), 36 * k, 2);
    lab.forEach((o, j) => {
      if (levels[j] < 0) return;
      card.append("text").attr("x", x(o.mz)).attr("y", up(o.y) - 7 * k - levels[j] * 11 * k)
        .attr("text-anchor", "middle").attr("font-size", 9.5 * k)
        .attr("fill", pc).attr("font-weight", 600)
        .text(`${o.ppm >= 0 ? "+" : ""}${o.ppm.toFixed(0)} ppm`)
        .on("mousemove", e => showTip(tt, e, `<b>slot error</b><br>${o.ppm >= 0 ? "+" : ""}${o.ppm.toFixed(1)} ppm`, 150))
        .on("mouseleave", () => hideTip(tt));
    });

    // theoretical comb down: every stick carries its M+k label under its tip
    // (height differences stagger them); when slots crowd, alternate a second
    // row rather than dropping labels
    const t = p.theoretical, pitch = (r - l) / Math.max(t.x.length - 1, 1);
    const row2 = pitch < 26 * k ? j => (j % 2) * 11 * k : () => 0;
    t.x.forEach((mz, j) => {
      card.append("line").attr("x1", x(mz)).attr("x2", x(mz))
        .attr("y1", mid).attr("y2", down(t.y[j]))
        .attr("stroke", pc).attr("stroke-width", 2.2)
        .on("mousemove", e => showTip(tt, e,
          `<b>M+${t.offsets[j]}</b><br>${xlabel}: ${mz.toFixed(4)}<br>rel.: ${t.y[j].toFixed(1)}%`, 180))
        .on("mouseleave", () => hideTip(tt));
      card.append("text").attr("x", x(mz))
        .attr("y", Math.min(down(t.y[j]) + 11 * k + row2(j), bot - 3 * k))
        .attr("text-anchor", "middle").attr("font-size", 9 * k)
        .attr("fill", C.gMid).text(`M+${t.offsets[j]}`);
    });

    const yax = card.append("g").attr("transform", `translate(${l},0)`)
      .call(d3.axisLeft(y).tickValues([-100, -50, 0, 50, 100]));
    yax.select(".domain").attr("stroke-width", C.swAxis);
    yax.selectAll("text").attr("font-size", 10 * k).attr("fill", C.gDark);
    if (col === 0)
      card.append("text").attr("transform", "rotate(-90)")
        .attr("x", -(top + bot) / 2).attr("y", l - 44 * k)
        .attr("text-anchor", "middle").attr("font-size", 11 * k)
        .attr("fill", C.gMid).text(ylabel);
    const xax = card.append("g").attr("transform", `translate(0,${bot})`)
      .call(d3.axisBottom(x).ticks(5).tickFormat(d3.format(".2f")));
    xax.select(".domain").attr("stroke-width", C.swAxis);
    xax.selectAll("text").attr("font-size", 10 * k).attr("fill", C.gDark);
    card.append("text").attr("x", (l + r) / 2).attr("y", bot + 32 * k)
      .attr("text-anchor", "middle").attr("font-size", 11 * k)
      .attr("fill", C.gMid).text(xlabel);
  });

  window.__chartReady = true;
})();
