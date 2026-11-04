# Shared visual language

For template maintenance use templates/tokens.css and templates/micro.js.
Prefer existing colors, typography, tooltips and formatters for consistency.

Choose a scale by its meaning: sequential for ordered magnitude, diverging for
a meaningful center, categorical for classes. State the scale and legends;
missing observations remain visually distinct from measured low/zero values.
Check contrast and distinguish categories beyond color alone.

Make quantities easiest to compare through position or length where the view
allows it. Extra size/color channels should answer another useful question,
with a scale that does not exaggerate the variation.

Use units on axes and sensible precision for the measurement. Keep typography,
labels and spacing readable at the intended output size; journal requirements
come from the destination's current instructions, not a universal font floor
or remembered page width.

Keep canvas marks, SVG overlays and hit testing in one coordinate system.
Use exported values and legends that agree with the data contract. Record
source, transforms, thresholds and meaningful counts with the figure.
