# VarPep software architecture (updated downstream-only selection)

VarPep uses a staged Python/SQLite design to keep **complete evidence catalogs** independent from optional peptide candidate selection.

1. **Reusable indexes.** `build-index` creates a SQLite representation of transcripts and exons from GTF and transcript/protein FASTA. `build-rediportal-index` and `build-snp151-index` create separate lookup indexes.
2. **Evidence collection.** `call` reads a name-sorted BAM using `pysam`, aggregates variants and co-observed local events from informative sequencing fragments, and merges worker results into persistent SQLite state. Neither `--exclude-known` nor `--nonsynonymous-only` removes events during this stage.
3. **Complete catalogs.** The pipeline retains the supported variant and local transcript-haplotype catalogs, including events that do not alter the encoded protein. Annotation adds database evidence and ordered per-event interpretations.
4. **Optional peptide selection.** `--exclude-known` and `--nonsynonymous-only` restrict **only** the events/haplotypes contributing to peptide generation. The first uses strand-compatible REDIportal substitutions and SNP151 position-matched SNVs; the second selects events that individually change a translated protein. Both may be combined.
5. **Translation and peptide windows.** Selected haplotypes are projected onto annotated transcript coding sequences; altered protein sequences yield candidate MHC-I and MHC-II peptide windows.
6. **Processing predictions.** Pepsickle is applied to protein sequences associated with retained MHC-I peptide candidates. These predictions do not establish HLA binding or antigen presentation.

The resulting local phasing is limited to the span of informative RNA-seq fragments; it does not determine full-length transcript haplotypes. Annotations are not proof of germline/somatic origin or tumor-specific editing. Refer to the source for exact SQLite schema and event eligibility conditions.
