#!/usr/bin/env bash
set -euo pipefail
# Edit reference paths before running. Use a NEW out-prefix for each run.
python -u varpep.py call \
  --bam SAMPLE.name.bam --name-sorted \
  --reference GRCh38.primary_assembly.genome.fa \
  --transcripts gencode.pc_transcripts.fa \
  --proteins gencode.pc_translations.fa \
  --index gencode_varpep.sqlite \
  --rediportal-index rediportal_v3_hg38.sqlite \
  --snp151-index snp151_hg38.sqlite \
  --out-prefix SAMPLE_varpep_filtered \
  --include-indels --threads 8 \
  --exclude-known --nonsynonymous-only
# Variants and haplotypes are still exported in full.
