/* chart-forge micro library — the scaffolding every figure used to hand-roll.
 * Attaches window.forge; injected by build.py before the chart script.
 * Style source of truth is templates/tokens.css — forge.C reads it at load time,
 * so figures stay in sync with tokens without touching chart code.
 * d3 is referenced lazily (inside functions) so --bundle none / 3dmol flows
 * can still use forge.C and forge.journal(). */
window.forge = (function () {
  "use strict";
  const css = getComputedStyle(document.documentElement);
  const v = n => css.getPropertyValue(n).trim();
  const px = pt => Math.round(pt * 96 / 72 * 10) / 10;

  const C = {
    font: v("--font"),
    ink: v("--ink"), gDark: v("--g-dark"), gMid: v("--g-mid"),
    gSoft: v("--g-soft"), gFaint: v("--g-faint"), gGhost: v("--g-ghost"),
    grid: v("--c-grid"), threshold: v("--c-threshold"),
    accent: v("--c-accent"), unmod: v("--c-unmod"), link: v("--c-link"),
    xref: v("--c-xref"), gap: v("--c-gap"),
    ns: v("--c-ns"), pOnly: v("--c-p-only"), fcOnly: v("--c-fc-only"),
    bIon: v("--c-b-ion"), yIon: v("--c-y-ion"), aIon: v("--c-a-ion"),
    cIon: v("--c-c-ion"), zIon: v("--c-z-ion"), xIon: v("--c-x-ion"),
    fontMono: v("--font-mono"),
    cat: [v("--cat-1"), v("--cat-2"), v("--cat-3"), v("--cat-4"),
          v("--cat-5"), v("--cat-6"), v("--cat-7")],
    swAxis: parseFloat(v("--stroke-axis")) || 0.75,
    swGrid: parseFloat(v("--stroke-grid")) || 0.5,
  };

  /* ---- tooltip card (tokens.css --tip-*) ---- */
  function tip() {
    return d3.select("body").append("div")
      .style("position", "fixed").style("opacity", 0).style("pointer-events", "none")
      .style("background", v("--tip-bg")).style("color", v("--tip-fg"))
      .style("padding", v("--tip-pad")).style("border-radius", v("--tip-radius"))
      .style("font-size", v("--tip-fs")).style("font-family", C.font)
      .style("line-height", "1.5").style("z-index", "10");
  }
  function showTip(t, e, html, w) {
    t.html(html).style("opacity", 1)
      .style("left", Math.min(e.pageX + 14, window.innerWidth - (w ?? 270)) + "px")
      .style("top", e.pageY - 10 + "px");
  }
  const hideTip = t => t.style("opacity", 0);

  /* ---- journal preset resolution (null on screen builds → keep token defaults).
   * Resolved once and cached; axes/grid/colourbar consult it automatically, so
   * any chart using forge primitives is journal-compliant for free. ---- */
  let _J;   /* undefined = unresolved, null = no preset */
  function journalNow() {
    if (_J !== undefined) return _J;
    const J = window.JOURNAL;
    if (!J) { _J = null; return _J; }
    _J = {
      name: J.name,
      mmSingle: J.mmSingle, mmDouble: J.mmDouble,
      W: Math.round(J.mmSingle * 96 / 25.4),
      Wdouble: Math.round(J.mmDouble * 96 / 25.4),
      fsMin: px(J.ptMin), fsTick: px(J.ptTick), fsLabel: px(J.ptLabel),
      gridlines: !!J.gridlines, axisBox: !!J.axisBox, ticksIn: !!J.ticksIn,
      font: J.font || C.font, panelLabels: J.panelLabels || "",
    };
    return _J;
  }

  /* figure geometry: resolves width from the journal preset when present,
   * else the screen default. Returns {W, H, k, J} — scale margins by k. */
  function fig(o) {
    const J = journalNow();
    const W = J ? (o && o.double ? J.Wdouble : J.W) : (o && o.W) || 1054;
    const k = W / 1054;
    const H = Math.round((o && o.H) * k);
    return { W, H, k, J };
  }

  /* ---- gridlines behind data (call before mark layers, or re-render per zoom).
   * Journal presets that ban gridlines (nature) suppress them here. ---- */
  function gridV(g, x, o) {   // vertical lines: needs plot height h and top t
    if (journalNow() && !journalNow().gridlines) return;
    const { h, t, n = 8 } = o;
    x.ticks(n).forEach(tick => g.append("line")
      .attr("x1", x(tick)).attr("x2", x(tick)).attr("y1", t).attr("y2", t + h)
      .attr("stroke", C.grid).attr("stroke-width", C.swGrid));
  }
  function gridH(g, y, o) {   // horizontal lines: needs plot width w and left l
    if (journalNow() && !journalNow().gridlines) return;
    const { w, l, n = 6 } = o;
    y.ticks(n).forEach(tick => g.append("line")
      .attr("x1", l).attr("x2", l + w).attr("y1", y(tick)).attr("y2", y(tick))
      .attr("stroke", C.grid).attr("stroke-width", C.swGrid));
  }

  /* ---- axes with token defaults; fonts and tick direction follow the journal
   * preset automatically (ACS ticks point inward). Labels optional. ---- */
  const axFs = () => journalNow() ? journalNow().fsTick : 11;
  function axBottom(g, x, o) {   // {y, n=8, tickSize, label, ly}
    const J = journalNow();
    const sel = g.append("g").attr("transform", `translate(0,${o.y})`)
      .call(d3.axisBottom(x).ticks(o.n ?? 8)
        .tickSize(o.tickSize ?? (J && J.ticksIn ? -4 : 3)));
    sel.select(".domain").attr("stroke-width", C.swAxis);
    sel.selectAll("text").attr("font-size", axFs());
    if (o.label) sel.append("text")
      .attr("x", (x.range()[0] + x.range()[1]) / 2).attr("y", o.ly ?? 34)
      .attr("font-size", J ? J.fsLabel : 12)
      .attr("fill", C.gMid).attr("text-anchor", "middle").text(o.label);
    return sel;
  }
  function axLeft(g, y, o) {     // {l, n=6, tickSize, label, lx}
    const J = journalNow();
    const sel = g.append("g").attr("transform", `translate(${o.l},0)`)
      .call(d3.axisLeft(y).ticks(o.n ?? 6)
        .tickSize(o.tickSize ?? (J && J.ticksIn ? -4 : 3)));
    sel.select(".domain").attr("stroke-width", C.swAxis);
    sel.selectAll("text").attr("font-size", axFs());
    if (o.label) sel.append("text").attr("transform", "rotate(-90)")
      .attr("x", -(y.range()[0] + y.range()[1]) / 2).attr("y", o.lx ?? -40)
      .attr("font-size", J ? J.fsLabel : 12)
      .attr("fill", C.gMid).attr("text-anchor", "middle").text(o.label);
    return sel;
  }

  /* ---- bottom-centre colourbar: ramp gradient + lo/hi/label (tokens.md: render
   * the ACTUAL scale, never a swatch approximation). Pass colour stops directly so
   * sequential and diverging ramps share one code path. ---- */
  function cb(svg, colorFn, o) {   // {lo, hi, W, y=null, w=220, title, fmt=(f)=>f.toFixed(1), n=21}
    const J = journalNow();
    const fs = J ? J.fsMin : 10;
    const w = o.w ?? 220, cbX = o.cx ?? (o.W - w) / 2, cbY = o.y;
    const defs = svg.append("defs");
    const grad = defs.append("linearGradient").attr("id", "cb-" + Math.random().toString(36).slice(2, 7));
    const n = o.n ?? 21, f = o.fmt ?? (x => x.toFixed(1));
    d3.range(n).forEach(i => grad.append("stop")
      .attr("offset", i * (100 / (n - 1)) + "%")
      .attr("stop-color", colorFn(o.lo + (o.hi - o.lo) * i / (n - 1))));
    svg.append("rect").attr("x", cbX).attr("y", cbY).attr("width", w).attr("height", 10)
      .attr("fill", `url(#${grad.attr("id")})`).attr("rx", 2);
    svg.append("text").attr("x", cbX - 8).attr("y", cbY + 9).attr("text-anchor", "end")
      .attr("font-size", fs).attr("fill", C.gMid).text(f(o.lo));
    svg.append("text").attr("x", cbX + w + 8).attr("y", cbY + 9)
      .attr("font-size", fs).attr("fill", C.gMid).text(f(o.hi));
    if (o.title) svg.append("text").attr("x", cbX + w / 2).attr("y", cbY + 24)
      .attr("text-anchor", "middle").attr("font-size", fs).attr("fill", C.gMid).text(o.title);
    return { cbX, cbY, w };
  }

  /* ---- journal panel label (Nature convention: lowercase bold a, b, c) ---- */
  function panel(svg, ch, x, y, fs) {
    return svg.append("text").attr("x", x).attr("y", y)
      .attr("font-weight", 700).attr("font-size", fs ?? 11).attr("fill", C.ink).text(ch);
  }

  /* ---- journal preset resolution (null on screen builds → keep token defaults).
   * Width comes from mm @96dpi; fonts from the preset's pt minimums. ---- */
  /* public journal(): journalNow() with explicit column choice —
   * journal({double:true}) → double-column width in W (the preset's mm). */
  function journal(o) {
    const J = journalNow();
    return J && o && o.double ? { ...J, W: J.Wdouble } : J;
  }

  /* ---- greedy label dodge: items must be pre-sorted by x; returns a level
   * (0..levels-1) per item, or -1 when every level is taken within minGap.
   * Shared by peak-label stacks, volcano callouts, anything column-adjacent. ---- */
  function dodgeX(xs, minGap, levels = 2) {
    const placed = Array.from({ length: levels }, () => []);
    return xs.map(x => {
      const lv = placed.findIndex(col => col.every(px => Math.abs(px - x) > minGap));
      if (lv < 0) return -1;
      placed[lv].push(x);
      return lv;
    });
  }

  /* ---- 1-D brush zoom controller (feature-map contract): drag-select to zoom,
   * tiny selections ignored, selection cleared after fire, reset restores.
   * The caller's onZoom([a,b]) re-renders and owns window.__zoomState. ---- */
  function zoomBrush(svg, o) {   // {extent, x, onZoom, minPx=8}
    const brush = d3.brushX().extent(o.extent)
      .on("end", ev => {
        if (!ev.selection) return;
        const [a, b] = ev.selection;
        const el = svg.select(".brush-zoom");
        if (b - a < (o.minPx ?? 8)) { el.call(brush.move, null); return; }
        o.onZoom([o.x.invert(a), o.x.invert(b)]);
        el.call(brush.move, null);
      });
    // beneath the plot content: the brush overlay has pointer-events:all and
    // would otherwise swallow every mark's real-mouse hover. Redraws append
    // marks back on top, so the layering stays correct across zooms.
    svg.insert("g", ":first-child").attr("class", "brush-zoom").call(brush);
    return {
      reset(fullDomain) { o.x.domain(fullDomain); o.onZoom(fullDomain); },
    };
  }

  /* ---- panelStack: vertical cascade of panels sharing one x scale.
   * Each panel gets its own <g>; draw(g, top, ctx) returns the panel height
   * (fixed-height panels just return the number). `static: true` panels draw
   * once and keep their marks across stack.redraw() — the sequence title
   * does not follow zoom; every non-static panel is wiped and re-drawn.
   * Atoms take absolute tops (prsmTitle/errorPanel), so draw receives the
   * panel's top directly and the <g> stays untransformed. ---- */
  function panelStack(svg, panels, o) {
    const gap = (o && o.gap) ?? 16;
    const ctx = { W: (o && o.W) || 1054, k: (o && o.k) || 1 };
    const gs = panels.map(() => svg.append("g"));
    const heights = panels.map(() => 0);
    let top = (o && o.top) || 0;

    function drawPanel(i) {
      const g = gs[i];
      g.selectAll("*").remove();
      const h = panels[i].draw(g, top, ctx);
      heights[i] = typeof h === "number" ? h : 0;
      top += heights[i] + (panels[i].gap ?? gap);
    }
    panels.forEach((_, i) => drawPanel(i));

    return {
      /* total height after the last panel's trailing gap — usable as the
       * svg height's bottom pad */
      get H() { return top; },
      redraw() {
        top = (o && o.top) || 0;
        panels.forEach((p, i) => {
          if (p.static) { top += heights[i] + (p.gap ?? gap); return; }
          drawPanel(i);
        });
      },
    };
  }

  /* ---- legend row: swatch circles + labels, horizontal pitch or vertical stack ---- */
  function legend(svg, entries, o) {   // entries [{label, color}], o {x, y, pitch=110, dy=18, fs=11, vertical}
    const J = journalNow();
    const fs = o.fs ?? (J ? J.fsMin : 11);
    const g = svg.append("g");
    entries.forEach((e, i) => {
      const ex = o.x + (o.vertical ? 0 : i * (o.pitch ?? 110));
      const ey = o.y + (o.vertical ? i * (o.dy ?? 18) : 0);
      g.append("circle").attr("cx", ex).attr("cy", ey).attr("r", 5).attr("fill", e.color);
      g.append("text").attr("x", ex + 10).attr("y", ey + 3.5)
        .attr("font-size", fs).attr("fill", C.gDark).text(e.label);
    });
    return g;
  }


  /* ---- number formatting presets (d3.format, resolved lazily) ---- */
  function fmt(kind) {
    const m = { si: ".3s", int: ",", f2: ".2f", ppm: "+.2f" };
    return d3.format(m[kind] || kind);
  }

  return { C, tip, showTip, hideTip, gridV, gridH, axBottom, axLeft, cb, panel, fig, journal, dodgeX, zoomBrush, legend, panelStack, fmt };
})();
