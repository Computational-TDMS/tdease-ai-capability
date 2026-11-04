# Ion conventions and localization

Choose series from the acquisition and supported kernel conventions.

| Activation | Typical series |
|---|---|
| CID / HCD | B,Y |
| ETD / ECD | C,Z |
| Mixed or other activation | Inspect metadata and specify supported series |

The kernel exposes plain Z for ETD/ECD; ZDot and ZPrime are rejected. Do not
substitute radical-ion mass conventions from another program. Use the kernel
for theoretical masses instead of reconstructing its atomic-offset table.

When max charge is omitted, the kernel uses declared precursor charge, with
fallback/source recorded in maxChargeSource. Check this metadata when it is
missing or doubtful; choose overrides from acquisition evidence. A larger
theoretical search space also increases opportunities for accidental matches.

For a chain of length n, a modification at position p affects N-terminal ions
numbered at least p and C-terminal ions numbered at least n-p+1. Observed ions
that distinguish placements support localization; bracketing can come from one
series or complementary series. Report the unresolved span/set when observed
predictions remain equivalent.

Fragmentation and modification retention depend on activation, charge and
chemistry. The usual behavior of one method is a prior, not a guarantee for a
particular scan. UniMod records can supply per-site neutral-loss predictions:

    tdeasenexus-toolbox unimod-get --title <modification>

A neutral loss can support a group without localizing its residue.
