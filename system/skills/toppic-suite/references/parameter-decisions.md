# Parameter Decisions

Always query the deployed Node before writing exact flags. Use this reference to choose
semantics, not to override the Node contract.

## Units

- TopFD peak/deconvolution error is commonly expressed in `m/z`.
- TopPIC fragment/precursor matching tolerance is commonly expressed in `ppm`.
- PrSM cluster/proteoform error and unexpected mass-shift ranges are expressed in `Da`.

Reject values such as `10ppm` when the binary expects the bare integer `10`.

## TopFD choices

Set these from the acquisition method or cited workflow:

- activation method
- precursor isolation window
- maximum charge
- maximum monoisotopic mass
- MS1 and MS2 signal-to-noise ratios
- peak/deconvolution error
- missing-MS1 behavior

Never assume the acquisition has MS1 scans. Probe the converted mzML's scan
levels before choosing flags: an MS2-only method (no MS1 spectra) requires
`topfd -o` (missing-level-one), and then TopPIC has no TopFD feature file to
read, so it also requires `toppic -x` (no-topfd-feature). One data fact drives
both flags; record it in config.yaml rather than rediscovering it from a
TopFD failure.

Recover acquisition facts from metadata and supplied protocols. Use deployed
Node defaults for routine settings when they fit the objective, and record
them. Ask for missing experimental facts or a scientific choice that materially
changes the intended analysis, rather than requiring confirmation of every flag.

Do not assume binary defaults reproduce a cited paper. Obtain reported settings
from the user-supplied protocol or paper record and store them as run-specific
configuration; they are not global skill defaults.

## TopPIC choices

Choose these from the objective, supplied protocol, existing Run and deployed
Node contract. Reuse prior decisions; ask only about unresolved choices that
materially change the user's intended analysis:

- matching tolerance
- proteoform/PrSM cluster tolerance
- maximum number and range of unexpected shifts
- allowed proteoform types
- allowed N-terminal forms
- fixed or variable modifications
- EVALUE versus FDR cutoff
- target-only versus target-decoy search

FDR requires a valid decoy strategy. Do not label an EVALUE cutoff as FDR. If a paper
reports only one FDR threshold, inspect the matching TopPIC version/method details before
deciding whether it applies at spectrum, proteoform, or both levels; record the mapping.

Published TopPIC settings are run-specific. Record matching tolerance,
clustering, FDR, unexpected-shift range and maximum shifts in the run manifest
rather than embedding one dataset's values in the skill. When clarification is
needed, explain the relevant reference value and trade-off; the Wilson & Zhou
chain is an example protocol, not a mandatory parameter questionnaire.

## Literature reproduction

Record all compatibility changes separately from cited parameters:

```text
reported method -> deployed Node version -> required compatibility adjustment
```

Do not silently replace a literature setting with a newer binary default.
