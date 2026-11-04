# Reference knowledge

The first prior pass is defined in ../SKILL.md. This reference explains how to
use the results, rather than adding another lookup sequence.

| Source | Useful context |
|---|---|
| UniProt | Sequence, processing, entry annotations and known sites |
| dbPTM / iPTMnet | Protein/site precedent, enzymes and record provenance |
| ELM / cohort | Sequence motif and annotated context |
| UniRef | Family members and related protein context |
| BLAST | Homolog alignment and possible site transfer |

Use ptm-database for atomic or combined lookup. Registered UniProt tools also
support get, query, blast and map. Choose a full sequence for family context
or a local window for a site question; follow up when the initial results
leave a useful uncertainty.

Keep exact-protein annotations, homolog annotations and motif matches distinct.
For site transfer, retain taxonomy, aligned residues and coordinate mapping.
An empty or unavailable source is not evidence that the sample lacks a PTM.
Record source identity/date when supplied.

Domain memory and user suggestions can also propose or rank explanations.
Label uncited ideas as AI judgment or user hypothesis; do not invent a
literature citation. Reference knowledge does not measure this sample.

Reuse cached, source-matched records when available. A disconnected frontend
may still permit direct REST access; no network means use local mirrors/cache
and report remote sources as unavailable. BLAST needs network.
