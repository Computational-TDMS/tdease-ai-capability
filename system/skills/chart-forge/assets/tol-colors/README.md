# Paul Tol colour schemes (CVD-safe qualitative sets)

Source: github.com/Descanonge/tol_colors (pip `tol-colors` 2.2.0) — packaging of Paul
Tol's SRON technical note on colour schemes for data. Free to use; schemes designed for
colour-vision-deficiency safety. `tol_colors.py` is vendored verbatim as upstream
reference; the D3-consumable materialization lives at `templates/tol-colors.js`
(TOL_BRIGHT / TOL_HIGH_CONTRAST / TOL_MEDIUM_CONTRAST / TOL_VIBRANT / TOL_MUTED /
TOL_LIGHT).

Positioning vs our tokens: `tokens.css --cat-1..7` (Okabe-Ito) is the default
categorical set; reach for TOL_MUTED (10 colours) or TOL_VIBRANT when a figure needs
more/distinct categorical series. `light` is NOT CVD-safe — pale/print contexts only.
