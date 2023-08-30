# STATUS

Done:
- file names reflect desired conventions
- rebuilt fastq from mis-mapped `bam` files

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

```bash
sbatch --array=1-103%10 scripts/align_remade_fastq.sh
# 6436542
# to rerun:,

sbatch --array=10,66,69,70,71,72,73,74,75,76,77,9 scripts/align_remade_fastq.sh

sbatch --array=1-200%10 scripts/align_fastqs.sh
# 6436555

sbatch --array=0,65,66,73,76,80,98 scripts/align_remade_fastq.sh
```

```bash
DATADIR="${PWD}/INPUT/FASTQ_FILES_RENAMED/"
FASTQS=($(ls ${DATADIR}/*_R1_001.fastq.gz))
IDS=($(basename -a ${FASTQS[@]%_R1_001.fastq.gz}))
echo ${IDS[@]} | \
    tr ' ' '\n' | \
    grep -n -f <(cat <(ls INPUT/FASTQ_FILES_RENAMED/ | cut -d 'R' -f 1 | sed 's/_$//g' | sort -u ) | grep -v -f <(ls OUTPUT/UNMERGED_BAMS/ | cut -d '.' -f 1 )) | \
    sed 's/:/\t/g' | \
    awk '{print $1-1,$2}' | \
    cut -d ' ' -f 1 | \
    tr '\n' ','

sbatch --array=141 scripts/align_fastqs.sh
```

## Truncated file fixing
One pair of `fastq` files appeared to contain missing reads or improperly named pairs, such that
`bwa` would give an error and quit. I used `trimmomatic` to extract properly paired reads:

```bash
# set directory
cd /vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/FASTQ_FILES_RENAMED

# create renamed symlink
ln -s /vf/users/CARDPB/data/HBCC/ILM_WGS_data/FASTQ_FILES/H6-98455_S21_L003_R1_001.fastq.gz 82053_L003_R1_001.fastq.gz
ln -s /vf/users/CARDPB/data/HBCC/ILM_WGS_data/FASTQ_FILES/H6-98455_S21_L003_R2_001.fastq.gz 82053_L003_R2_001.fastq.gz

# capture all paired reads with trimmomatic
module load trimmomatic

java -jar ${TRIMMOJAR} PE -threads 2 \
    82053_L003_R1_001.fastq.gz \
    82053_L003_R2_001.fastq.gz \
    82053_L003_R1_paired.fastq.gz \
    82053_L003_R1_unpaired.fastq.gz \
    82053_L003_R2_paired.fastq.gz \
    82053_L003_R2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 MINLEN:36

# rename files to use only paired R1 and R2
rm 82053_L003_R1_001.fastq.gz
rm 82053_L003_R2_001.fastq.gz
rm 82053_L003_R1_unpaired.fastq.gz 
rm 82053_L003_R2_unpaired.fastq.gz 
mv 82053_L003_R1_paired.fastq.gz 82053_L003_R1_001.fastq.gz
mv 82053_L003_R2_paired.fastq.gz 82053_L003_R2_001.fastq.gz
```


```bash
# get-missing-FASTQ_FILES_RENAMED

REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
DATADIR="${PWD}/INPUT/FASTQ_FILES_RENAMED/"
FASTQS=($(ls ${DATADIR}/*_R1_001.fastq.gz))
IDS=($(basename -a ${FASTQS[@]%_R1_001.fastq.gz}))

N=-1
FILES=($(ls INPUT/FASTQ_FILES_RENAMED/*R1_001.fastq.gz))
for fn in ${FILES[@]%_R1_001.fastq.gz}; do
    let N=${N}+1
    fs=$(basename $fn)
    if [ ! -f "OUTPUT/UNMERGED_BAMS/${fs}.bam" ]; then
        echo "sbatch --array=${N} scripts/align_fastqs.sh # for ${fs}"
    fi
done
```


```bash
REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
DATADIR="${PWD}/INPUT/BAM_FILES"
BAMS=($(ls ${DATADIR}/*.final.bam))
IDS=($(basename -a ${BAMS[@]%.final.bam}))

FILES=($(ls INPUT/REMADE_FASTQ_FROM_BAM/*R1_001.fastq.gz))
for fn in ${FILES[@]%_R1_001.fastq.gz}; do
    fs=$(basename $fn)
    if [ ! -f "OUTPUT/REMADE_BAMS/${fs}.bam" ]; then
        N=$(grep -n $fs convert.txt | cut -d ':' -f 1)
        let N=${N}-1
       echo "sbatch --array=${N} scripts/align_remade_fastq.sh # for ${fs}"
    fi
done


sbatch --array=76 scripts/align_remade_fastq.sh # for 81948
sbatch --array=80 scripts/align_remade_fastq.sh # for 81949
sbatch --array=98 scripts/align_remade_fastq.sh # for 81962
sbatch --array=65 scripts/align_remade_fastq.sh # for 82020
sbatch --array=73 scripts/align_remade_fastq.sh # for 82047
sbatch --array=66 scripts/align_remade_fastq.sh # for 82050
```

# MERGE fastqs
```bash

IDS=(81926-FASTQ
81933
81937
81942
81946
81950
81955
81956
81957
81960
81970
81971
81972
81975
81976
81978
81979
81986
81988
81989
81993
81994
81995
81996
81998
81999
82000
82001
82006
82019
82022
82031
82034
82043
82044
82046
82053
82055
82058
82062
82064
82065
82067
82072
82076
82077
82083
82084
82085
82130)

sbatch scripts/merge_bams.sh ${IDS[@]:0:10}     # 7082928
sbatch scripts/merge_bams.sh ${IDS[@]:10:10}    # 7082931
sbatch scripts/merge_bams.sh ${IDS[@]:20:10}    # 7083033
sbatch scripts/merge_bams.sh ${IDS[@]:30:10}    # 7083035
sbatch scripts/merge_bams.sh ${IDS[@]:40:10}    # 7083039

sbatch scripts/merge_bams.sh 81994              # 7161139
```


# `fastq` file splitting
To split paired-end `fastq` files into 4 parts (for each R1 and R2):

```bash
module load seqkit
seqkit split2 \
-1 ${id}_R1_001.fastq.gz \
-2 ${id}_R2_001.fastq.gz \
-p 4 \
-e .gz \
-O ${id}_split 
id='81962'
```