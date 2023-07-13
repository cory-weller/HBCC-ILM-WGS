# STATUS

Done:
- file names reflect desired conventions
- rebuilt fastq from mis-mapped `bam` files
- generated tiny toy `fastq` files

To do:
- align original `fastq` to `hg38` reference
- merge `bam` across lanes
- align remade `fastq` (from `bam`) to `hg38` reference
- generate `gVCF`
- generate joint called `VCF`


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
See [`build_test_data.sh`](scripts/build_test_data.sh).
The script generates small subset `fastq` files with 100k rows each.
Test files located in `INPUT/TEST`.

    tinyBAMfastq_R1_001.fastq.gz
    tinyBAMfastq_R2_001.fastq.gz
    tinyfastq_L001_R1_001.fastq.gz
    tinyfastq_L001_R2_001.fastq.gz


## Snakemake workflow
