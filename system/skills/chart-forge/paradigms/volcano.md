# Volcano design

Template-maintenance reference. The shipped statistical report uses the shared
volcano panel; this is not a standalone template route. The statistical
producer owns tests and adjustment.

Label x as the stated log fold change and y as negative log10 of the actual
p-value or adjusted value. State contrast direction, adjustment, sample counts
and thresholds. Make axis clipping, capped significance or excluded values
explicit. Centering the effect axis helps comparison; choose ranges for the
data and purpose without hiding asymmetry.

Threshold membership can use categorical color, with a legend. Additional
channels need distinct meanings and readable scales; if size encodes amount,
make area proportional and supply a size legend. Omit channels that add little
information or exaggerate a narrow range.

Prioritize labels for the question and avoid overlap at the final size.
Tooltips and selection can expose detailed identity and estimates. Density
rendering may reduce overplotting, with the displayed aggregation disclosed.

Presence-only features may be outside the producing test's estimable set.
Report them separately according to the method, rather than inventing finite
fold changes or significance. A pattern or threshold crossing does not supply
causal interpretation or validate the experimental design.
