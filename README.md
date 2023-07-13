# HBCC-ILM-WGS

## Prepare input

Generated `.tsv` files for the `.bam` and `.fastq` inputs, along with
their desired new IDs: [`bam_samples.tsv`](bam_samples.tsv) and
[`fastq_samples.tsv`](fastq_samples.tsv)

Input data includes
- 154 total samples
- 104 bams to convert back to fastq then align
- 50 sammples across 400 FASTQ files (4 lanes each)


## Fix `bam` files

Most of the `bam` files complained about lack of `EOF` block when checked with `samtools quickcheck`. 
See [`bam_check.md`](REPORTS/bam_check.md) for a list of which files did and did not show missing `EOF` block.
The files were converted to `fastq` with [`bam_to_fastq.sh`](scripts/bam_to_fastq.sh).
Briefly, files were converted to uncompressed bam using  `samtools view -u`.
Then read pairs extracted with `samtools fastq`. 
This also converted file IDs from `OLD_ID` to `NEW_ID` as defined in [`bam_samples.tsv`](bam_samples.tsv).



```bash
sbatch --array=0-103%20 scripts/bam_to_fastq.sh
# jobid 4385862
```

Four jobs had stalled and were re-ran to completion:
```bash
sbatch --array=73,81,84,85 scripts/bam_to_fastq.sh
# jobid 4417955
```


## File renaming

The goal: rename `fastq` files with new IDs while Keeping the original files intact.
To not use disk space for redundant files, I instead created symlinks replacing `OLD_ID` with `NEW_ID`.
This was done with [`create_symlinks.py`](scripts/create_symlinks.py) and only for the `fastq` files
(as the `bam` files were renamed by [`bam_to_fastq.sh`](scripts/bam_to_fastq.sh)).

```bash
python3 scripts/create_symlinks.py
```


## Tiny test data
Generating small subset fastq files with 100k rows
```bash
# parameters
nrow='100000'
testdir='INPUT/TEST'
mkdir -p ${testdir}

# Tiny version of (original) fastq
infilestem='INPUT/FASTQ_FILES_RENAMED/82084_L001'
outfilestem="${testdir}/tinyfastq_L001"

zcat "${infilestem}_R1_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R1_001.fastq.gz"

zcat "${infilestem}_R2_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R2_001.fastq.gz"

# Tiny version of regenerated (from bam) fastq
infilestem='INPUT/REMADE_FASTQ_FROM_BAM/81925'
outfilestem="${testdir}/tinyBAMfastq"

zcat "${infilestem}_R1_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R1_001.fastq.gz"

zcat "${infilestem}_R2_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R2_001.fastq.gz"
```

## Snakemake workflow
