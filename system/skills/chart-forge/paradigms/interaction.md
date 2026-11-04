# Interaction for template maintenance

Use interaction when it helps inspect scientific values or relationships.
Routine analysis uses the shipped adapter; adding interactions is not an
analysis prerequisite.

Keep scale meanings and series identities stable during hover, selection and
zoom. Tooltips reveal exact values and source IDs; essential axes, units and
deciding values remain readable without hover. Selection can emphasize related
marks without obscuring the original encoding.

Choose controls for the view: zoom for dense regions, series toggles for
overlap, or linked selection for related panels. Avoid adding every control to
every chart. Zoomed views need a visible range and reset; labels adapt to the
visible domain.

For spectra, link peaks, assignments and sequence coverage using verified
fragment-to-residue mapping. For structures, preserve residue coordinates
between sequence and 3D selections. For mirror panels, synchronize the axis.

Production assets use shipped D3 or 3Dmol. Vendored viewers offer design
examples, not runtime dependencies. Set window.__chartReady when rendering is
complete for the browser audit; expose additional state when export or
inspection needs it.

Check useful interactions and the static view at the intended size. A figure
must remain interpretable when printed or when a tooltip is unavailable.
