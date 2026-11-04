# Adapting data for a figure

Read the selected template's data.json and build.json. Match its actual
contract; illustrative values demonstrate shape, not expected results.

The domain producer owns parsing, scientific transformations and provenance.
Preserve source identity, coordinate conventions, units, selection/exclusions
and relevant uncertainty. Use an existing source reader rather than rebuilding
mzML or TopPIC parsing in a plotting script.

A field's missing-value semantics come from its source. Represent unavailable
plotted measurements as null and distinguish them from measured zero. If a
matrix uses zero as a non-detection sentinel, handle that explicitly in the
adapter and record the transform.

Do not impose tidy rows, field naming or aggregation rules that contradict
the selected template. Dense interval views may need exact source rows;
distribution views may use documented bins or sampling. Record counts and
transformations so display reduction is visible.

Supply meta.source and the relevant transform, thresholds, units and sample
counts. Keep the plotted contract alongside the output; the figure should
remain traceable without recovering hidden browser state.

If no shipped view represents the question, provide inspectable results and
identify the missing view. A new reusable template is maintenance, not a
required detour during analysis.
