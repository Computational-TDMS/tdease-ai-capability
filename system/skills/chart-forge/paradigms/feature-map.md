# Feature and family maps

Template-maintenance reference. Runtime routes are examples/feature-map/ and
examples/ms1-family-map/; read their data.json for the actual contract.

Preserve the observation's geometry. A recorded RT interval can be a horizontal
bar with an apex marker; an apex-only measurement is a point. The shipped
feature view uses RT in minutes and average m/z. A family view represents
spacing/member relations, not the same object.

Choose color, size or facets for the question, with explicit legends.
Intensity and charge are possible variables, not universal assignments.
Draw important marks legibly without silently omitting weaker observations.

Use SVG or canvas according to measured density and interaction performance,
rather than a fixed mark-count law. Canvas views need axes/labels and usable
hit testing; PNG export preserves their marks. Linked panels and zoom share
data coordinates and provide a reset.

Show filtering, row counts, transforms and source in metadata. An RT shift,
an unusual charge band or a family spacing invites investigation; the map
alone does not identify a chromatography fault or chemical increment.
