# Mass-shift tool recipes

Invocation details only. ../SKILL.md owns investigation choices and the early
reference context. Use CLI help for unfamiliar options.

## Refine an uncertain MS1 frame

Reuse an existing source-matched refinement. When mass/isotope registration
is an open question and suitable MS1 observations are available:

    tdeasenexus-toolbox refine-shift <source.mzML> --proforma <declared-ProForma> \
      --mz <precursor-m/z> --charge <z> --rt-min <min> --rt-max <max> \
      --mod-mass <unexpected-Da> --prior-shift <unexpected-Da>

Take source, m/z, charge, RT and declared terminal/modification frame from the
target itself. Supply the RT window explicitly; the default is near run start.
--prior-shift is in the post-declared basis; refinedShiftAfterDeclared is the
residual in that basis, while refinedMass is the total neutral mass.
--observed can reuse a previously extracted envelope JSON.

Here --mod-mass supplies the unassigned residual to the theoretical mass;
--prior-shift supplies the calibration prior, not candidate chemistry. The
numeric-residual path uses the declared baseline's isotope shape and cannot
establish the residual's composition. For a complete chemical candidate, supply
its complete ProForma without adding that residual again; compare identity
models with spectral-evidence under a fixed frame.

Inspect ambiguity, warnings, matchedPeaks, pearson and evaluations (coverage,
residuals, fittedDa and isotope offsets). Compare the fit with source/family
masses. A weak or competing anchor retains the raw spread and a provisional
frame. fittedDa is a shared
axis offset, not candidate evidence; fit precision is not systematic accuracy.

## Localize or decompose a residual

OpenSearch proposalInput.argv binds source ProForma, scan and sibling MSAlign
when available. missingInputs and limitation describe gaps, not absence of
spectra elsewhere. Use the source-linked argv or the explicit forms:

    tdeasenexus-toolbox resolve-shift <*_ms2.msalign> --scan <scan[,scan...]> \
      --proforma <source-ProForma> [--library builtin|none|<json>]

    tdeasenexus-toolbox resolve-shift <*_ms2.msalign> --scan <scan[,scan...]> \
      --sequence <sequence> --shift <Da> --mods <fixed-mods-json> \
      [--library builtin|none|<json>]

When the bracket shift and a reliable refined frame differ, use --sequence,
--shift and the fixed components to state that frame explicitly. Choose
informative source-linked scans; add independent scans to resolve uncertainty
or disagreement. Report the tested set and its limitations.

| Option | Meaning |
|---|---|
| --library omitted or builtin | Small common-PTM proposal set plus a free residual |
| --library none | Localize one free residual without a decomposition lens |
| --library <json> | Source-derived modificationLibrary or an explicit lens; entries {name,residues,location?,positions?} |
| --max-components | Default 3; adjust for the plausible model space |
| --no-free-shift | Diagnostic library-only closure test |
| --overshoot-window-da | Default 0.06 Da for near-closure proposals |
| --details | Supporting details for returned candidates, without changing the search |
| --ion-types | Normally inferred from activation; unknown activation needs explicit series |

Pool only precursor-compatible sibling scans under the same activation.
Read result.candidates, component intervals, mass error, freeMass, ties and
warnings. Keep source digests and the full output in artifacts.
closedByLibrary describes accounting, not identity. Generated proforma can be
reused directly, but a row with freeMass still has an unresolved component.

## Query reference priors

ptm-database owns query syntax, source coverage and returned records.
../SKILL.md owns how reference context informs the investigation. Reuse those
instructions and source-matched results rather than adding a second query
sequence here.

## Search mass candidates

    tdeasenexus-toolbox unimod-search --mass <residual-Da> --tolerance <Da> \
      [--residues <letters>] [--limit <N>]

    tdeasenexus-toolbox unimod-solve --target-mass <residual-Da> --tolerance <Da> \
      [--max-modifications <1..4>] [--allowed-site <R[:position]>] [--limit <N>]

Use the stated residual basis and measured uncertainty. Start with a focused
search; use combinations, alternative terminal frames or broader lenses when
they address the live explanation. Processing can create new terminal sites,
so constrain the processed model rather than only the original sequence.

For several selected targets, unimod-batch-solve accepts an --input JSON array
of {id, mass, tolerance} and the same combination/site filters as unimod-solve.
It loads the database once; use it when those queries are independent. Keep
each row's id tied to its proposed backbone/processing basis and uncertainty.
The batch computes combinations for supplied targets; it does not establish
which processing model is present or require enumerating every possible basis.

### Search with a conditional fixed component

Let T be observed neutral mass minus the chosen reference's exact mass.
If additional components F are pinned, the solver searches T - sum(F).
Components already included in the reference are excluded from T and must
not also be supplied as pins. Backbone/processing choices belong to that
reference; --fixed-mod pins database modification identities and counts.

    tdeasenexus-toolbox unimod-solve --target-mass <T> --tolerance <Da> \
      --fixed-mod <TITLE|UNIMOD:n|USER:n> --max-modifications <total-count>

Supply the target INCLUDING the pins; do not manually subtract their mass
and then pass them again. --max-modifications includes pinned occurrences.
Alternatively search the already-subtracted remainder without those pins,
then reassemble and check the whole model. Allowed sites constrain placement
of the whole combination; use the chosen backbone's coordinates and available
sites. A modification pin alone is not a measured site assignment.

Batch input rows may add fixedMods: [<reference>, ...] to their {id, mass,
tolerance}. Global --fixed-mod applies to every row; row pins add to it, and
repeated entries add counts. Keep branches requiring different site filters
in separate calls. Read result.queries[].solution for each batch result;
result.queries[].fixedMods records its effective pins.

A pin makes a conditional search, not evidence that the component is present.
Retain its origin and uncertainty, and compare credible unpinned alternatives
when the pin is unconfirmed. unconstrainedBaseline relaxes query filters but
keeps the pins; it does NOT test a model without that component. Rejection of
an unknown or filtered-out pin is an input error, not a chemical no-hit.

Keep stdout and stderr separate:

    tdeasenexus-toolbox unimod-solve ... > solve.json 2> solve.warn

Read result.candidates for solves; mass search uses result.massMatches
(not text-query result.matches), with title, deltaMono, signedErrorDa and
specificities. Inspect matchCount, returnedCount, truncated, warnings and
unconstrainedBaseline when provided. A parse failure or truncated head is not
a no-hit.

Save full results and declared filters. Rank a working set using source
constraints, priors and distinguishing predictions; group equivalent
predictions and identify the untested remainder. A database match, including a
unique mass closure, remains a proposal.

## Fragment matching and tool location

The kernel normally infers ion series from activation and the fragment charge
ceiling from precursor charge. Override with a stated physical/method reason.
Read resolved maxCharge and maxChargeSource when reporting match counts.
Raw mzML matching uses isotope envelopes; `spectral.py measure` supplies the
comparison as an evidence card, and `spectral.py show` the isotope windows.

Hand spectral-evidence the complete candidate set as a `models.json`
(`[{"id","label","proforma"}, ...]`, one bare sequence) and keep the card it
returns: the card is the measurement, and it states per-ion state, isotope
offset and the peaks two models both claim. Report a comparison from the card,
not from coverage or match counts alone.

In the product shell, tdeasenexus-toolbox is on PATH. Skill scripts locate the
packaged kernel automatically. Use installed entry points and Workspace outputs.
