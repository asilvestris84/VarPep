# VarPep revised downstream-only selection

- `--exclude-known`: only restricts peptide generation; dbSNP151 is position-based and REDIportal requires A>G/T>C compatible with catalog strand. No exclusion at the BAM evidence-collection stage.
- `--nonsynonymous-only`: only restricts peptide generation to haplotypes whose constituent events independently change the encoded protein. Joint-only effects of individually silent variants are not selected with this flag.
- The complete fragment-supported variant catalog is retained, including synonymous and noncoding SNVs and indels.
- The haplotype catalog includes all globally supported raw transcript haplotypes for validated transcript models, even those without protein consequences. Haplotype enumeration still follows the existing fragment and transcript rules.
- Variant columns: `nonsynonymous_in_any_transcript` and `potential_neoantigen_candidate` (yes/no/unresolved).
- Haplotype columns: `event_classifications`, `event_rsids_ordered`, and `contains_potential_neoantigen_event`; semicolon-separated annotations preserve event order.
- Pepsickle runs only on proteins associated with retained MHC-I peptide candidates. VarPep does not itself perform HLA binding prediction.

## Interpretation and limitations

- dbSNP151 position match is **not** proof of germline origin; dbSNP151 absence is **not** proof of a somatic mutation.
- REDIportal absence is **not** proof of tumor-specific RNA editing; editing requires strand-aware and ideally matched DNA/normal validation.
- A variant-level `nonsynonymous_in_any_transcript=yes` indicates at least one evaluated expressed transcript context; different isoforms can have different effects.
- `potential_neoantigen_candidate` is a heuristic annotation, not experimental antigen presentation.
- The two downstream flags act independently: combining both restricts peptides to events that meet both criteria.
- The workflow requires indexed REDIportal and SNP151 databases for the `--exclude-known` flag.

## Validation status

Syntax compilation and limited smoke tests were performed. Full RNA-seq integration, performance and peptide-processing tests were not possible in this environment because pysam and the required BAM/reference/index inputs are unavailable. Run an end-to-end validation on a small `--region` with a fresh `--out-prefix` before production.
