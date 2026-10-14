# Annotation and downstream-only filtering

## Complete catalogs

The output `PREFIX.variants.tsv` includes all retained variants and reports `nonsynonymous_in_any_transcript` and `potential_neoantigen_candidate` (yes/no/unresolved). The output `PREFIX.haplotypes.tsv` retains supported local transcript haplotypes and reports `event_classifications`, `event_rsids_ordered`, and `contains_potential_neoantigen_event`. Semicolon-separated event annotations follow the same order as the haplotype `events` field.

## REDIportal and dbSNP151

REDIportal provides strand-dependent known editing annotations. SNP151 provides positional SNV matches and rsIDs. **The current SNP151 matching is by genomic position, not exact REF/ALT.** A database match is evidence of prior annotation, not proof of germline origin; absence is not proof of a somatic event. REDIportal absence does not establish tumor specificity.

## Selection for peptide output only

`--exclude-known` restricts peptide generation by excluding events matching REDIportal or dbSNP151 under the implemented criteria. Both indexes are required. `--nonsynonymous-only` restricts peptide generation to events individually protein-changing in the evaluated transcript context. The options can be combined. Neither option removes entries from the exported variant or haplotype catalogs. Pepsickle processing is restricted to retained MHC-I candidate contexts.

A putative neoantigen annotation is **not** experimental neoantigen evidence. Validation requires matched DNA/normal tissue, allele-specific confirmation where possible, and downstream immunological studies.
