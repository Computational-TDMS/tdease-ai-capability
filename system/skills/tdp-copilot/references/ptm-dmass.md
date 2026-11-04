# Mass basis and processing

Read when mass accounting or isotope registration affects candidate selection.
Shared investigation principles come from the system prompt; domain routing
is in ../SKILL.md and commands are in tool-recipes.md.

## Establish the basis

Total delta is observed intact neutral monoisotopic mass minus the theoretical
mass of the stated backbone. A TopPIC unexpected shift is usually a residual
after terminal processing and named PTMs, not that total delta. Recover the
source's terminal form and assigned modifications before comparing models.

For a stated backbone, the mass balance is
`observed mass = reference mass + sum of signed component changes`.
The residual is what remains after the components already accounted for.
If the reference mass changes by `dB`, the same observation's residual changes
by `-dB`. Re-express both the target and its site constraints in the new basis;
a no-hit in one basis is not a no-hit in another.

A net composition change can contain negative element counts when atoms or
residues are removed. Nonnegative atom-count and standalone-molecule valence
filters test a standalone species, not a general signed difference between
species. Judge chemical feasibility on the complete proposed species.

Use parse-proforma for exact model masses and composition. Domain memory can
supply candidate ideas; attach their origin and distinguish recollection from
a queried reference or sample measurement.

## Treat accuracy as measured

Deconvolution may misregister an isotope; large proteins may have an inferred
rather than observed monoisotopic peak. Keep neutral Da, m/z and ppm distinct.
An intact-mass ppm window scales with intact mass, not the small residual.

refine-shift supplies a provisional source frame, not an external mass standard.
Its fit interval is statistical precision, not systematic accuracy. Compare it
with the source/family spread and preserve calibration disagreements. Weak or
missing refinement permits exploration with wider uncertainty.

## Processing and candidate space

### N-terminal assignments from a search

Treat a reported N-terminal form as part of the search's proposed model, with
support conditional on that search and its observations. Recover the actual
sequence interval, terminal assignment, fixed/variable PTMs and unexpected
shifts before interpreting a residual. Read the TopPIC output field contract
through system:tdp-opensearch, resource references/toppic-output.md; loading
the search-execution Skill is unnecessary.

When sequence context, reference knowledge or spectral conflicts motivate an
alternative terminal model, retain the original assignment and derive the
alternative from the same observed mass. Update the backbone mass, available
sites, coordinate mapping and residual together; count processing and assigned
components exactly once. A change of notation is not a new observation.

Equal intact-mass balances do not establish equal chemical identities or
fragment predictions. Use relevant terminal fragments or other discriminating
observations when available; missing terminal coverage may leave processing
and modification alternatives unresolved. A better numerical closure alone
does not prove a terminal state. No automatic enumeration of every terminal
form is required, and a no-hit alone does not show that TopPIC was wrong.

### Candidate selection

Backbone, processing and modification choices affect reference mass and
available sites. When source facts support competing models, state each basis
and the corresponding residual explicitly. Select queries and measurements
for their ability to distinguish those models.

A combination search mixes biological modifications, processing, sequence
changes and preparation reagents. Rank candidates in context rather than
discarding all combinations together. Save the full output, including filters
and truncation, and choose a working set by feasibility and distinguishing
predictions; exhaustive individual measurement is unnecessary.

Near-isobaric models may need MS1 composition or MS2 discrimination. Mass
agreement, even a unique database hit, is compatibility rather than identity.
Keep an unresolved signed residual when no defensible model accounts for it.
