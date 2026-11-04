# Figure build and verification

Runtime figures use the selected domain adapter or build.py; the engineering
contract defines package inputs and outputs. Scientific assets inline shipped
dependencies/data for offline use.

For template/compiler maintenance, edit the package, build a representative
figure and inspect it at the intended size. Exercise relevant interactions
and exports. Fix observed clipping, unreadable labels, wrong encodings or
runtime errors, then verify the affected behavior again.

scripts/audit.py --browser is the full maintenance gate.
scripts/check-figures.py checks produced files for self-containment; it does
not replace visual inspection or scientific validation.

Pure SVG views have in-page SVG export. Canvas/hybrid/WebGL need PNG through:

    node scripts/export.mjs --in <figure.html> --out <export-stem>

Keep HTML, exact plotted data and provenance together. Requested static
exports retain the scientific meaning and readable labels.
