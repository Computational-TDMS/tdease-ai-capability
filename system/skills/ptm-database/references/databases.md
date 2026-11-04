# Reference source scope

Choose a source for the question it can answer. The initial prior pass is an
investigation decision; this catalog adds no lookup order. scripts/query.py
owns access and scripts/fetch.py owns mirror provisioning.

| Source | Useful question | Interpretation limit |
|---|---|---|
| UniProt | Sequence, processing and annotated sites on the named entry | Check isoform, evidence and coordinates |
| dbPTM mirror | Site/type precedents and sequence-window cohorts | Coverage depends on the downloaded types and mirror date |
| iPTMnet mirror | PTM sites, linked enzymes and source scores | Source scores are not sample identification probabilities |
| ELM mirror | Matches to motif-class patterns | Regex matches omit website structural/compartment filters |
| UniRef / BLAST | Related sequences for follow-up | Homolog annotations require alignment before site transfer |
| PhosphoSitePlus | Additional site/regulatory annotations | User-supplied licensed data; not bundled or redistributed |
| Proteoform Atlas | Prior intact-proteoform observations | Match source sequence/context, not just mass |

Use the report's source identity, version/date and actual lane status.
Distinguish a matching record, an empty query, an unavailable source and a
skipped query. No annotation is not evidence of biological absence.

The BLAST adapter returns hit identities and summary scores, not alignment
coordinates. Retrieve an alignment before transferring a site or interpreting
an unaligned region. In flow, homolog annotation summaries cover only hits with
local mirror records; the report also retains the retrieved hit/member lists.

The local mirror normalizes source type spellings through its vocabulary.
dbPTM windows are 21 residues centered at index 10; iPTMnet counts deduplicated
sites, not papers. Query the mirror rather than recreating a manual download
or grep cascade during analysis.

External resource pointers:
- UniProt: https://www.uniprot.org/
- PhosphoSitePlus licensing: https://www.phosphosite.org/staticLicensing
- Proteoform Atlas: http://repository.topdownproteomics.org/

Access terms and remote interfaces can change. Check the current source when
obtaining new data; old catalog ratings, record counts and successful test
queries are not guarantees of coverage or access.
