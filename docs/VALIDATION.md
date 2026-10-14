# Validation checklist

**Status:** Python syntax compilation and limited smoke tests only; updated end-to-end calling and Pepsickle processing on a real BAM remain unvalidated.

1. Build all three indexes using the same reference assembly.
2. Run a small `--region` with **four distinct output prefixes**: no optional filters; `--exclude-known`; `--nonsynonymous-only`; both.
3. Compare `variants.tsv` and `haplotypes.tsv` across all runs. Counts, identities and database annotations should be identical under otherwise identical parameters.
4. Check that `peptides.tsv` and Pepsickle outputs are appropriately restricted by each option.
5. Check a negative-strand T>C event annotated by REDIportal and a known SNP151 rsID; verify per-event labels for two-event haplotypes.
6. Confirm peptide sequence correspondence to altered proteins and the absence of duplicated counting due to run-state reuse.
7. Record runtime, peak memory, input BAM characteristics and reference versions. These checks are not substitutes for an independent accuracy benchmark.
