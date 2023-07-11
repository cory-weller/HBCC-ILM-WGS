# HBCC-ILM-WGS

## Prepare input

Generated `.tsv` files for the `.bam` and `.fastq` inputs, along with
their desired new IDs: [`bam_samples.tsv`](bam_samples.tsv) and
[`fastq_samples.tsv`](fastq_samples.tsv)

Input data includes
- 154 total samples
- 104 bams to convert back to fastq then align
- 50 sammples across 400 FASTQ files (4 lanes each)

## File renaming
For simplicity of processing via snakemake, rename files to include the
desired CARD_ID contained in the above sample `stsv` files.

Keeping the original files intact (and not using disk space for redundant
files), I instead created symlinks replacing `OLD_ID` with `NEW_ID`.

See [`create_symlinks.py`](scripts/create_symlinks.py)

## Snakemake workflow

