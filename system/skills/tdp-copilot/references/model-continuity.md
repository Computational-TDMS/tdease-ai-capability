# Continuing and revising models

Reuse current source-linked artifacts and the user's corrections. Carry forward
the source/scan IDs, coordinate frame, terminal processing, assigned components,
candidate origins, relevant measurements and unresolved conflicts.

Check the whole model when changing one component. Source assignments,
fragment-supported localization and mass-only proposals have different status.
Revising a supported component calls for an explanation tied to the evidence;
a new mass closure alone does not invalidate its localizing fragments.

Keep a bounded working set including the user's requested models. Store broad
solver results in artifacts. Group observationally equivalent alternatives,
test representative predictions, and state what was grouped or left untested.
Equivalence is relative to the stated observations and settings. Preserve
distinct model identities when grouping; a different observation may separate
them and reopen the comparison.
An unresolved interval stays an interval until a measurement separates sites;
enumerate placements when that particular comparison needs them.

Use generated ProForma where available and check coordinates when constructing
or changing models; the measurement adapter validates supplied ProForma inputs.
A complete composition may still have unresolved sites; an open
residual is not a complete chemical model. Preserve complete decomposition
paths rather than silently mixing their components.

Report precursor incompatibility separately from an explicitly requested MS2
comparison when the tool accepts that complete model. Continuing the question
does not require repeating discovery, the prior pass or unchanged measurements.
