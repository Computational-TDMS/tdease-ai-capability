# Measurement choices

Select by the uncertainty. These methods are alternatives and complements,
not stages that every investigation must complete.

## Complete-model MS2 comparison

Use complete ProFormas including processing, fixed tags and co-modifications.
Compare on the same source, preprocessing, ion series, charge range, tolerance
and calibration treatment. The raw mzML route uses isotope-envelope matching;
TopFD msalign contains deconvoluted neutral masses and supports a different,
single-model annotation route.

`spectral.py measure` validates and measures the supplied models in one pass, and
returns the verdict as a card. A direct kernel call is useful when only the
match tables are needed. Sum scans only when precursor identity and activation
justify pooling; preserve scan-level conflicts that pooling could conceal.

Read the card's `compositionDeltaDa` first: it says how far apart the models
actually are. A near-isobar pair (differences of a few Da, or a mass defect
close to an integer number of neutrons) has offset isotope combs, so one
model's monoisotopic fragment can sit on the other's M+2 — check
`contestedPeaks` before counting either model's matches as independent
support. Aggregate superiority can be informative without establishing which
chemistry or site caused it.

## MS1 composition

Extract the physical feature once and fit complete theoretical compositions
against that shared observation. Site placements with identical composition
have the same MS1 prediction. Use shared calibration/registration constraints;
candidate-specific flexibility must not manufacture a better match.

For a family, compare source-linked members using the same method and each
member's complete composition. Inspect fit quality and stability across
members. Recurrence concerns the measured family, not transfer of identity or
localization to unmeasured members.

## Registration and calibration

Integer isotope assignment and continuous mass-axis calibration are different
uncertainties. Use source information to constrain them and apply comparable
treatment to alternatives. If registration is ambiguous, compare the plausible
registrations explicitly rather than choosing a favorable offset per model.

A systematic signed-error pattern can motivate checking calibration, isotope
assignment, ion convention or interference. Inspect charge states and relevant
independent observations before attributing it to one cause.

## MS1 and MS2 envelope scores are not one scale

The MS2 fragment matcher (apex-anchored monotone walk, missing slots counted
as zero, Pearson over all envelope slots, correlation floor) and the MS1
feature evaluator (nearest peak per theoretical slot, Pearson over matched
pairs only, padded-zero explicitly rejected) are different methods with
deliberately different definitions of "how well the envelope fits". Both are
valid for their own question; their correlation numbers are not comparable
evidence scores, and a threshold calibrated on one does not transfer to the
other. Unifying the two algorithms underneath would be a methods change
requiring its own validation — until then, report which method produced a
number whenever a correlation is cited.

## Breakpoints, losses and diagnostic ions

A processing breakpoint needs observed fragments that distinguish the proposed
boundary; missing fragments alone are inconclusive. Complementary-pair sums can
check consistency with the precursor under the actual ion convention.

For a loss or diagnostic ion, specify predicted m/z, charge and tolerance, and
check precursor compatibility, activation and resolving power. Observation can
support a chemical group; non-observation depends on sensitivity and the scan
set. Use UniMod records or grounded composition to derive predictions.

## Inspecting isotope windows

Use `show` for the selected fragment envelopes of one MS2 scan. Underneath,
`envelope-regions` returns all raw peaks within each theoretical window with
isotope-slot and ppm annotations: extraction, not a new correlation measure. A
window can look well populated and still be rejected on the full-slot Pearson,
so read the kernel's verdict and the window's slot count together.

`show --ms1` addresses whole-precursor envelopes, not individual MS2 ions.
Adapters render recorded measurements; figures add no independent evidence.
