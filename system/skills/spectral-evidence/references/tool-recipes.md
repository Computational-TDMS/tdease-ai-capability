# Spectral measurement commands

One adapter entry, three verbs. Run it from the Workspace and put outputs in
the question's analysis directory.

    SPECTRAL="$TDEASE_SYSTEM_ASSETS/skills/spectral-evidence/scripts/spectral.py"

## measure — the criterion

    uv run --script "$SPECTRAL" measure --mzml <file.mzML> \
      --scan <scan[,scan...]> --models <models.json> --out <card.json>

Scans select by comma list, `--scans-input` (inline JSON or a file, for
large batches), an RT window (`--rt-min/--rt-max`) or a scan range
(`--scan-min/--scan-max`); `--sum-scans` pools compatible scans
explicitly. Every model is measured on one shared extraction, whatever
the scan set.

`models.json` is an array of `{"id", "label", "proforma"}` — every model
complete (bare sequence, terminal forms, all modifications) and sharing one
bare sequence. One model gives a compatibility audit; two or more give the
comparison. The kernel runs one envelope-matching pass per scan for every
model on the same extracted peaks; the full match tables go to
`<card>-evidence.json` (megabytes — do not read them by hand), and the verdict
comes back as the card, which is both the stdout and `--out`. The default
output is also the full figure set: the overview figure plus the separating
isotope windows (the selection `show` would draw), written beside the card —
`--no-figures` keeps the measurement only. Nothing downstream needs
rendering by hand; named ions go through `show`.

Stdout is the card, so nothing has to be assembled from sidecars:

- `observation` — the file, scans, ion series, activation, charge ceiling and
  tolerance actually used, and the criterion.
- `models[]` — each model's own numbers (`matchedIons`, `coverage`,
  `discriminatingMatched`, `medianAbsErrorPpm`, `p95AbsErrorPpm`,
  `uniqueCleavages`), plus its proforma so a later `show` needs only the card.
- `verdict.compositionDeltaDa` — the mass distance between the models, from
  the kernel's parse of each complete model.
- `verdict.contestedPeaks` — observed peaks matched by more than one
  model/ion pair, with each claim's `isotopeOffset` and ppm.
- `verdict.ions[]` — the ions whose prediction differs, most separating first:
  `role`, `matchedBy`, `cells` (`cellsByScan` for a multi-scan run), `state`
  per model, `reportsSites`, `coversDiffSites`.
- `siteWindows` — the kernel's cleavage-window counts at the differing sites.
- `uninformative`, `warnings`, `next`, `audit`.

Options: `--ion-types` (otherwise the activation decides: ETD/ECD → C,Z;
HCD/CID → B,Y), `--activation`, `--max-charge` (otherwise the precursor's
declared charge — the physical ceiling), `--tolerance-ppm` (default 10),
`--focus-site`, `--min-rel-intensity`, `--ms-bin`.

## measure --ms1 — precursor envelopes

    uv run --script "$SPECTRAL" measure --ms1 --mzml <file.mzML> \
      --feature <LABEL> <MASS> <CHARGE> <RT_MIN> <RT_MAX> \
      --models <models.json> --out <card.json>

Repeat `--feature` for a family. The kernel aggregates the observed precursor
envelope over each RT window in one mzML traversal, and each model's
theoretical isotope envelope is generated at that charge. Same card shape,
different verdict: `verdict.features[].ranking` gives each model's composition
mass (`anchorMass`), its distance from the observed feature mass
(`anchorOffsetDa`), shape agreement (`pearsonR`) and apex delta. MS1 separates
compositions, not sites — two models differing only in *where* a modification
sits have the same precursor envelope. Use `--isotopes` to fix the isotope
count; 0 lets the kernel derive it.

## annotate — msalign leads

    uv run --script "$SPECTRAL" annotate --msalign <sample_ms2.msalign> \
      --scan <SCANS-value> --proforma <complete-model> --out <annotation.html>

Writes the figure, the fragment-match sidecar and a `*-card.json`. Preserve the
linked source scan. `--annotate 'mz=LABEL'` / `'mz1..mz2=LABEL'` adds attention
pointers — annotations, not measured data.

## show — the windows

    uv run --script "$SPECTRAL" show <card.json> [--ions <ion:charge,...>] [--ppm 10]

One panel per model per ion from `envelope-regions`, written beside the figure
as `<stem>.regions<i>.json`. Defaults to the card's separating ions plus two
controls; at most 12 windows per call. The charge is required per ion — a bare
ion expands to every charge state and stops being one window.

    uv run --script "$SPECTRAL" show --ms1 --refine <refine-shift.json> \
      [--refine <second.json>] --out <envelope.html>
    uv run --script "$SPECTRAL" show --ms1 --observed <feature.json> \
      --theory <a.json> --theory <b.json> [--fit <a-fit.json>] --out <overlay.html>

Renders an MS1 measurement that already exists. `--refine` draws the fitted
precursor envelope with its isotope-offset evaluations; `--observed` overlays
each complete model's theoretical envelope on the observed feature. Fit metrics
are bound to the pair by the fit's own recorded inputs — a proven mismatch is
rejected, missing provenance is reported as unverified.

---

# Appendix — kernel calls this skill does not wrap

Run `tdeasenexus-toolbox <command> --help` for the full contract. These return
kernel envelopes and are used directly.

| Call | Use |
|---|---|
| `fragment <msalign> --mode match --scan S --proforma M --evidence-out F` | single-model match when you need the tables yourself |
| `resolve-shift <msalign> --scan S ...` | shift discovery / localization (composition flow) |
| `envelope <mzML> --mass <M> --charge <z> --rt-min <a> --rt-max <b>` | extract one observed precursor envelope |
| `theory --proforma M --charge z` | theoretical isotope distribution |
| `fit-envelope --observed O --theoretical T` | Pearson fit of one pair |
| `refine-shift --observed O --mzml F ...` | refine a precursor mass from the observed envelope |
| `envelope-regions <mzML> --scan S --proforma M --ions <ion:charge>` | one isotope window's raw peaks and theoretical sticks (no verdict) |
| `spectrum <mzML> --mode scans --scans S[,S..]` | one scan's own peak list |
| `spectrum <mzML> --mode product --mz <m/z> --tolerance <ppm> [--scans S]` | census of one m/z across MS2 scans; **pass `--scans` to keep it about the scans you mean** |
| `select <mzML> ...` | scan selection |

`envelope-regions` is pure extraction: it returns every peak in the window with
its nearest theoretical slot, and passes no judgment. Pair it with a `measure`
card, never with a self-written matcher.

## Adjudicating a characteristic-loss or diagnostic-ion hypothesis

`measure` compares complete models and never annotates loss series. When a
candidate chemistry predicts a characteristic loss — an acyl chain leaving a
fragment, a labile group leaving the precursor, a diagnostic low-m/z ion —
test that position directly with the atomic extraction→theory→fit chain.
This is the move that breaks **mass degeneracy**: two candidate
modifications can agree at precursor mass to well inside any tolerance
(sulfation +79.9568 vs phosphorylation +79.9663, 0.0095 Da apart) while
predicting *different* characteristic losses (SO3 79.9568 vs H3PO4 97.9769,
18 Da apart) — the precursor position cannot separate them, the loss
positions can. Fit both members of the degenerate pair and the asymmetry
decides the chemistry.

The surviving species after the loss is what you theorize: a whole chain
minus the labile group (anchor at precursor neutral − loss), or a fragment
minus its acyl/mod (a b-ion anchors at peptide mono − 18.010565 − loss):

1. **Theory for the loss species.** Write the surviving fragment as a
   proforma carrying only what remains after the loss, and anchor at the
   *ion* mass — a b-ion is the peptide minus its C-terminal water, so
   `--anchor <peptide mono − 18.010565>`, then `--charge 1`:

       tdeasenexus-toolbox theory --proforma "<fragment+surviving mods>" \
         --charge 1 --anchor <ion-neutral-mass> > theory.json

   The composition only shapes the isotope distribution; the anchor places
   M0 where the ion is predicted.
2. **Observed peaks at that position.** Pull the MS2 scan's peaks and keep
   the ±1.6 Da neighbourhood of the prediction (enough for M0..M+2 at z=1):

       tdeasenexus-toolbox spectrum --mode window --ms-level 2 \
         --mz <predicted> --scan-min S --scan-max S <mzML>

   The window mode returns the scan's truncated peak list; filter to the
   neighbourhood and save `[{"mz":…, "intensity":…}]`.
3. **Fit.**

       tdeasenexus-toolbox fit-envelope --observed obs.json \
         --theoretical theory.json --match-window-ppm 30

Read the fit like any envelope: matched ≥ 2 isotope slots with a positive
Pearson is a measured species; a single lone peak has no envelope (Pearson
`null`) and supports nothing by itself. Fit the competing loss hypothesis the
same way and report the pair — the asymmetry is the adjudication. State the
prediction you tested, the slots matched, and the residuals; the fit output
is the artifact to quote.
