---
name: spectral-evidence
description: Measure MS1/MS2 spectral support for proteoform models, site placements and compositions.
---

# Spectral Evidence

Choose a measurement that could distinguish the stated models or answer the QC
question. Return observed differences and their limits.

One entry point, three verbs. Everything else in this skill is either a module
those verbs import or a kernel call listed in the appendix.

    uv run --script "$TDEASE_SYSTEM_ASSETS/skills/spectral-evidence/scripts/spectral.py" <verb> ...

| Verb | Question it answers |
|---|---|
| `measure --mzml F --scan S --models M.json` | Do these models actually differ, and on which ions? (card + overview + windows by default) |
| `measure --mzml F --ms1 --feature ... --models M.json` | Do their precursor envelopes differ — composition |
| `annotate --msalign F --scan S --proforma M` | Where are the leads in the deconvoluted spectrum? |
| `show CARD.json` | What do the separating windows actually look like? |

## The criterion, and the one thing not to do by hand

**mzML envelope matching is the criterion.** `measure` on the raw centroided mzML
runs the kernel's isotope-envelope matcher: the observed cluster is anchored on
its most abundant isotope, walked with a monotone cursor, missing slots count
as zero, and the full-slot Pearson must clear the kernel's correlation floor.
A hit on M+2 is a hit — report the isotope slot, correlation and ppm as
observed; do not invent a stricter monoisotopic-only standard.

**Do not decide existence yourself.** A monoisotopic-peak window written over a
peak list produces false negatives that contradict the criterion the rest of
the analysis uses. When a match surprises you, read its `state` in the card
(`missingEnvelope`, `missingApex`, `correlationRejected`, `unreasonableCharge`)
and look at the window with `show`. The card carries the kernel's measured
correlation even on a rejection, so there is nothing left to re-measure.

**msalign is the annotation view, not the criterion.** TopFD has already
deconvoluted the spectra, so a miss there carries no existence verdict —
collapse, deconvolution failure, filtering, interference and insufficient
signal all look the same. The accurate statement is that this deconvoluted
result provides no matching entry; check the raw data before concluding
anything from an msalign miss. Use `annotate` for leads; `resolve-shift`'s
sequence-graph search runs in this deconvoluted space.

## Reading the card

`measure` prints the card — the same bytes it writes to `--out`. Read it in this
order:

1. **`verdict.compositionDeltaDa`** — how far apart the models' complete
   chemistry is, from the kernel's parse of each full model. Two models ~2 Da
   apart are near-isobar: their isotope combs are offset, so one model's
   monoisotopic fragment sits where the other's M+2 sits.
2. **`verdict.contestedPeaks`** — observed peaks matched by more than one
   model/ion pair. These are not independent support for either model, and no
   single-model view can show them; a peak explained two ways is a normal
   competition, marked for inspection, not an elimination.
3. **`verdict.differences`** — predicted difference (`predictedDiffer`),
   identical prediction (`predictedIdentical` — this observation cannot
   separate the models there) and measured difference (`measuredDiffer`),
   stated separately. Resolvability depends on the measurement conditions the
   card reports; it is not pre-judged.
4. **`verdict.ions`** — the ions whose prediction differs, most separating
   first: ranked by how many models matched versus not, then match quality
   (correlation, error) and peak strength — never by the isotope slot a match
   landed in. Each row carries `reportsSites` (the residues whose
   modification state moves this fragment — an N-terminal c-ion of length n
   covers 1..n, a C-terminal z-ion covers the last n) and `coversDiffSites`.
   An ion whose span misses the differing site cannot be evidence about it,
   whatever the kernel says about its match.
5. **`uninformative`** — the remainder, counted: `identicalPrediction` (ions
   every model predicts identically; shared support only) and `noModelMatched`
   (no accepted envelope match anywhere — absence of acceptance, not evidence
   of absence unless observability is argued). The matched subset of
   identicalPrediction is surfaced as **`verdict.degenerateMatched`**: measured
   ions every model predicts at the same m/z, listed with
   `degenerateAmong` (the models they cannot separate). They are real
   measurements that ground span/window statements — "the whole differing
   block lies inside this ion's span" — and never a distinction between the
   models listed in `degenerateAmong`. `verdict.roleGlossary` states each
   row role's claim semantics in place.
6. **`next`** — capabilities the evidence leaves open, with the inputs each
   needs. It reports what this observation could not do; choosing what to do
   next belongs to the investigation.

Numbers you quote come from card cells, never from your own arithmetic. Every
matched cell carries the full tuple — `theoreticalMz`, `observedMz`,
`errorPpm`, `correlation`, `intensity`, `isotopeOffset` — in one m/z
convention. If you sum residue and modification masses yourself you will
eventually drop a proton or a charge; a number you computed that disagrees
with its card cell by anything is wrong, and the fix is to re-quote the cell,
not to re-derive it.

Coverage and match counts alone do not resolve sites. A single model's card
audits compatibility; it separates nothing.

## Looking at the windows

`measure` writes the default windows by default — the separating ions plus two
unmatched ones, one panel per model — alongside the overview figure, so one
measurement call is already the full figure set and nothing is left to render
by hand. `show` is for named ions:

    uv run --script .../spectral.py show <card.json> [--ions <ion:charge,...>]

One panel per model per ion: that model's theoretical isotope comb against the
raw peaks in the same window. A populated window in one model and an empty one
in the other is then visible rather than asserted. Each panel's subtitle
carries the kernel's verdict *and* the window's own slot count — they are
different statistics and reading either alone misleads (a window can look well
populated and still be rejected on the full-slot Pearson). Defaults to the
card's separating ions plus two windows with no accepted match, so absence is
visible too.

`show --ms1` renders an MS1 measurement that already exists (a refine-shift
result, or an observed feature against complete-model theories).

## Presenting the evidence in a reply

When a measurement grounds a claim, present the evidence in its working
order — the same order the investigation itself followed:

1. **Orientation — the overview figure** (`-figure.html`: annotated spectrum,
   per-model error maps, ppm-error histograms). It shows where across the
   sequence the models match and where they part, and thereby points at the
   decisive ions. It orients; it never carries the verdict alone — its
   histograms are descriptive, and a verdict stated from the error map alone
   has skipped the criterion.
2. **Solid evidence — the deciding ions' MS2 envelope windows** (the
   `-windows.html` of `measure`/`show` on the ions the orientation pointed
   to). The envelope match is the criterion; the populated-vs-empty panels
   of the same ion across models are the proof a reader checks.

Embed both with `![[path]]`, in that order — **the HTML files themselves**:
they are self-contained and render live in the reply (interactive,
auto-fitted). Never rasterize a figure to PNG for a reply: the export
freezes the interactive figure, spends a rendering round, and tempts panel
cherry-picking — if a figure must shrink to fit a raster, the ion selection
has already drifted. PNG/SVG exports are for external use, on explicit
request. Custom plots may supplement; they never replace either. A bare
`[[link]]` to a card file is a pointer, not a presentation — a reader who
cannot see the error map and then the deciding windows in the reply has not
been shown the evidence.

## Choosing by the question you have

`annotate → measure → show` is one available route, not a mandatory sequence.
With explicit candidate models and raw spectra, go straight to `measure`; when
the measurement answers the question, no extra routing round is owed — the
standard figures the measurement produced are still part of the answer and
follow the presentation rule below.

- **annotate** interprets the deconvoluted spectrum's matches and leads. It
  is not an existence criterion for anything.
- **measure** returns a comparable measurement between candidates under stated
  conditions. The name does not imply a completed statistical hypothesis
  measure.
- **show** expands existing observations and models, one panel per model per
  ion. It has no matching standard of its own — every verdict in a subtitle
  is the one `measure` already measured.

When the models themselves are not yet formed, or none of the current ones
explains the observation, that is hypothesis work — model construction and
reference knowledge, not another measurement.

Deep envelope collection exists to check that the fragment matches you are
relying on are solid, and to test a fragment's existence quickly — which
means it is collected on the SAME proteins and scans those matches came from:
envelope evidence gathered on other objects verifies nothing about them.

Batch scans with a comma list on `--scan`, or one call per scan. Select scans
for discriminating predictions and acquisition quality, using source-linked
records; there is no scan quota. A comparison spanning several proteins must
include each protein's OWN strongest spectrum — rank scans per protein from
the PrSM inventory, not within a single species. Additional scans are useful
when they address uncertainty or conflict.

## Questions and limits

| Question | Route |
|---|---|
| Do these models differ, and on which ions? | `measure --mzml --scan --models` |
| Do their compositions differ? | `measure --mzml --ms1 --feature ... --models` |
| What does the spectrum look like for this model? | `annotate --msalign` |
| What does the separating window look like? | `show <card>` |
| An existing MS1 fit or composition overlay | `show --ms1` |
| Family increment | `measure --ms1` per member, or a shared observed feature |
| Loss, diagnostic ion, one raw scan | `spectrum --mode scans` (appendix) |

Keep source IDs and coordinate frames explicit. Missing site detail need not
block composition comparison; unresolved chemistry may block a complete-model
calculation. Continue useful measurements whose inputs are ready.

`measure` requires every model complete (bare sequence, terminal forms, all
modifications) and sharing one bare sequence — an off-by-one backbone silently
invalidates every per-ion comparison. That shared-sequence constraint is a
capability limit, not evidence against a different-backbone hypothesis:
preserve the hypothesis and choose an applicable measurement rather than
altering its sequence to fit the interface. The figure adapter does not display
C-terminal modification tags; measure those models and report the missing view.

When a comparison establishes a mass difference but not its chemistry, the
measurement and its limits inform the ongoing investigation when chemical
interpretation is in scope. The unresolved question determines the next action;
completion of this measurement alone does not close that question. The
`--models` JSON and the card are the handoff — the card is the measurement.

## References on demand

- references/tool-recipes.md: the verbs in full, the kernel-call appendix, and the
  characteristic-loss/diagnostic-ion adjudication chain (theory --anchor → MS2
  window → fit-envelope) that breaks mass degeneracy at the loss positions.
- references/measurement-methods.md: measurement choices and comparability.
- references/evidence-bounds.md: interpretation and ambiguous evidence.
- references/fragment-ions.md: activation, ion conventions and localization.
- references/visualization.md: figure routes and available panels.
