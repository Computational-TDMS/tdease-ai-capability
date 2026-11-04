// Structure + modification sites workbench — the one chart-forge family
// rendered by 3Dmol (D3 cannot do molecules; see paradigms/structure-3d.md).
// Switches, all client-side over the embedded PDB:
//   colour mode   grey cartoon (unmodified = grey, tokens.md) <-> pLDDT
//                 confidence bins (AlphaFold standard, read from B-factors)
//   site toggles  per-site checkbox shows/hides sphere + cartoon tint + label
//   labels        numbered site labels on/off
//   spin / reset  idle rotation and fit-to-model
// Legend chips zoom to their site; clicking any residue inspects it.
(function () {
  "use strict";
  const { pdb } = window.DATA.data;
  const meta = window.DATA.meta;
  const sites = meta.sites;
  const C = forge.C;
  const siteColor = i => C.cat[i % C.cat.length];
  const plddtBins = meta.plddt || [];   // pre-binned in gen-data from CA B-factors

  const el = document.getElementById("chart");
  const narrow = el.clientWidth < 600;
  el.style.height = narrow ? "460px" : "620px";
  el.style.position = "relative";   // 3Dmol's canvas is position:absolute — keep it inside #chart, not over the page
  const viewer = $3Dmol.createViewer(el, { backgroundColor: "0xffffff" });
  const model = viewer.addModel(pdb, "pdb");
  const defaultZoom = narrow ? 1.1 : 1.35;

  /* ---------- state + single render path */
  const state = { mode: "grey", labels: !narrow, spin: false, visible: sites.map(() => true) };

  function applyView() {
    viewer.setStyle({}, { cartoon: { color: C.unmod, arrows: true }, stick: { hidden: true } });
    if (state.mode === "plddt") {
      plddtBins.forEach(bin =>
        viewer.setStyle({ resi: bin.resis }, { cartoon: { color: bin.color, arrows: true } }));
    }
    sites.forEach((s, i) => {
      if (!state.visible[i]) return;
      viewer.addStyle({ resi: s.resi }, {
        sphere: { radius: 2.6, color: siteColor(i) },
        cartoon: { color: siteColor(i) },
      });
    });
    viewer.removeAllLabels();
    if (state.labels) {
      sites.forEach((s, i) => {
        if (!state.visible[i]) return;
        const atoms = model.selectedAtoms({ resi: s.resi });
        const ca = atoms.find(a => a.atom === "CA") || atoms[0];
        if (ca) viewer.addLabel(`${i + 1} · ${s.label}`, {
          position: { x: ca.x, y: ca.y, z: ca.z },
          fontSize: 11, fontColor: C.ink,
          backgroundColor: "rgba(255,255,255,0.88)",
          borderColor: siteColor(i), borderThickness: 1.4, padding: 3,
          showBackground: true,
        });
      });
    }
    viewer.render();
  }
  applyView();
  viewer.zoomTo();
  viewer.render();
  viewer.zoom(defaultZoom, 600);

  /* ---------- toolbar (tokens-styled, above the canvas) */
  const btn = `cursor:pointer;font:600 11px ${C.font};padding:2px 10px;` +
    `border:1px solid ${C.gSoft};background:transparent;color:${C.gMid};border-radius:3px;`;
  const bar = document.createElement("div");
  bar.style.cssText = `display:flex;align-items:center;flex-wrap:wrap;gap:12px;` +
    `font:11px/1.6 ${C.font};color:${C.gDark};padding:6px 2px;` +
    `border-bottom:1px solid ${C.gGhost};margin-bottom:6px;`;
  bar.innerHTML =
    `<span style="color:${C.gMid}">colour</span>` +
    `<span style="display:inline-flex;gap:0">` +
    `<button data-mode="grey" style="${btn}border-radius:3px 0 0 3px;">grey</button>` +
    `<button data-mode="plddt" style="${btn}border-radius:0 3px 3px 0;">pLDDT</button></span>` +
    `<label style="cursor:pointer"><input type="checkbox" id="ptm-labels"${state.labels ? " checked" : ""}> labels</label>` +
    `<label style="cursor:pointer"><input type="checkbox" id="ptm-spin"> spin</label>` +
    `<button id="ptm-reset" style="${btn}">reset view</button>` +
    `<span style="color:${C.gMid}">sites</span>` +
    sites.map((s, i) =>
      `<label data-chip="${s.resi}" style="cursor:pointer;display:inline-flex;align-items:center;gap:4px;` +
      `border-bottom:2px solid ${siteColor(i)};padding-bottom:1px;">` +
      `<input type="checkbox" data-vis="${i}" checked>` +
      `<b style="color:${C.ink};cursor:zoom-in;">${i + 1}</b> ${s.label}</label>`).join("") +
    `<span id="ptm-inspect" style="margin-left:auto;color:${C.gDark}"></span>` +
    `<span id="ptm-plddt-key" style="display:none;width:100%;color:${C.gMid};">` +
    plddtBins.map(bin => `<b style="color:${bin.color}">■</b> ${bin.label}`).join(" · ") + "</span>";
  el.before(bar);

  /* colour-mode segment: active styling + re-render */
  const modeBtns = [...bar.querySelectorAll("[data-mode]")];
  function paintMode() {
    modeBtns.forEach(b => {
      const on = b.dataset.mode === state.mode;
      b.style.background = on ? C.gDark : "transparent";
      b.style.color = on ? "#fff" : C.gMid;
    });
    bar.querySelector("#ptm-plddt-key").style.display = state.mode === "plddt" ? "inline" : "none";
  }
  modeBtns.forEach(b => b.addEventListener("click", () => { state.mode = b.dataset.mode; paintMode(); applyView(); }));
  paintMode();

  bar.querySelector("#ptm-labels").addEventListener("change", e => { state.labels = e.target.checked; applyView(); });
  bar.querySelector("#ptm-spin").addEventListener("change", e => {
    state.spin = e.target.checked;
    try { viewer.spin(state.spin ? "y" : false); } catch { viewer.spin(!state.spin); }
  });
  bar.querySelector("#ptm-reset").addEventListener("click", () => {
    viewer.zoomTo();
    viewer.zoom(defaultZoom, 400);
    viewer.render();
  });

  /* site chips: checkbox = visibility, number = zoom to site */
  bar.querySelectorAll("[data-vis]").forEach(cb => cb.addEventListener("change", () => {
    state.visible[+cb.dataset.vis] = cb.checked;
    applyView();
  }));
  bar.querySelectorAll("[data-chip]").forEach(chip => {
    chip.querySelector("b").addEventListener("click", () => viewer.zoomTo({ resi: +chip.dataset.chip }));
  });

  /* click any residue → inspector readout */
  viewer.setClickable({}, true, (atom) => {
    if (!atom) return;
    const out = document.getElementById("ptm-inspect");
    if (out) out.innerHTML = `inspecting <b>${atom.resn}${atom.resi}</b> (chain ${atom.chain})`;
  });
  el.ondblclick = () => { viewer.zoomTo(); viewer.zoom(defaultZoom, 400); viewer.render(); };

  /* readiness: 3Dmol renders synchronously after render(); small paint delay */
  setTimeout(() => { window.__chartReady = true; }, 500);
})();
