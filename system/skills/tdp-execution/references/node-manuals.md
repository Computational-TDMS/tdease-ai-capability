# Node-specific constraints

Installed rules.smk and the shipped upstream manual own flags and outputs.
Use these notes for known interface traps, checking the deployed version.

## Informed-Proteomics

PbfGen creates a spectrum store, ProMex finds MS1 features, and MSPathFinderT
searches a FASTA under an explicit modification model. The chain can serve a
targeted TopPIC cross-check, but prior TopPIC execution is not an intrinsic
requirement of these tools.

Choose FASTA and modifications from the objective. A focused accession set and
small modification model reduce search cost while restricting the estimable
space. The bundled targeted library is an example chemistry set, not a
universal sample model. Ground compositions in a source or reference library;
untested chemistry remains a hypothesis regardless of search agreement.

### Output names for the documented 1.1.8305 interface

    sample.mzML → sample.pbf
                → sample.ms1ft
                → sample_IcTarget.tsv / sample_IcDecoy.tsv / sample_IcTda.tsv
                → sample.mzid

Outputs follow the input stem; -o is an output directory. ProMex uses
sample.ms1ft, not sample_ms1.ms1ft; mzIdentML is sample.mzid, not
sample_IcTda.mzid. Declare the outputs supported by the actual invocation.

PbfGen may reuse an existing .pbf; MSPathFinderT may reuse existing results
unless overwrite is enabled. Verify source/settings compatibility before
reuse. ProMex features can be reused with unchanged relevant parameters.

| Flag | Meaning / trap |
|---|---|
| -d | FASTA; -f instead means fragment tolerance |
| -feature | Feature file; -m is not its alias |
| -tda | Integer mode; 1 requests target and decoy, not the string on |
| -mod | Modification file; an empty value selects unmodified-only search |
| -overwrite | Controls reuse of existing search results |
| -t / -f | Precursor/fragment tolerances in ppm |
| -act | Activation enum; check the deployed binary's accepted values |
| -MinLength / -MaxLength | Sequence bounds |
| -MinMass / -MaxMass | Neutral-mass bounds in Da |

Modification records use composition,residues,opt|fix,position,name under
NumMods=N. Example: C2H2O,K,opt,any,Acetyl. Suite parameter files may use
DynamicMod syntax instead; read the selected file's format rather than mixing
the two. Search-space size changes runtime; no sample count guarantees success.

Read QValue and other actual output columns with their error-control scope.
Targeted agreement tests that search model, not independent chemical identity.

## TopFD / TopPIC

See toppic-suite. Important traps:
- Tolerances need the numeric format and units expected by the binary.
- Activation enums differ across tools.
- MS2-only acquisition requires paired TopFD -o / TopPIC -x handling.
- Preserve the sample_ms2.msalign → sample_ms2_toppic_* table relation.

## msconvert

Vendor raw to mzML conversion uses the deployed Node's centroiding, precision
and compression settings. The raw pipeline template includes msconvert;
existing converted inputs can bypass it. Retain the mzML.source.json sidecar
and check activation/scan levels when choosing downstream flags.
