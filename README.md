
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

FILES=($(ls INPUT/BAM_FILES_RENAMED/*final.bam))
for fn in ${FILES[@]%.final.bam}; do
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
id='81962'

module load seqkit
seqkit split2 \
-1 ${id}_R1_001.fastq.gz \
-2 ${id}_R2_001.fastq.gz \
-p 4 \
-e .gz \
-O ${id}_split 
```

```bash
sbatch scripts/splitbam.sh 81948
sbatch scripts/splitbam.sh 81962
sbatch scripts/splitbam.sh 82047
sbatch scripts/splitbam.sh 82050
```

```bash
sbatch scripts/align_split_fastq.sh 81948
sbatch scripts/align_split_fastq.sh 81962
sbatch scripts/align_split_fastq.sh 82047
sbatch scripts/align_split_fastq.sh 82050
```

parallel -j 1 sbatch scripts/align_split_fastq.sh :::  81948 81962 82047 82050 ::: 1 2 3 4

7455136
7455138
7455139
7455140
7455141
7455142
7455143
7455144
7455145
7455146
7455147
7455148
7455149
7455150
7455151
7455152

sbatch scripts/merge_bams.sh 81926-FASTQ
sbatch scripts/merge_bams.sh 81948
sbatch scripts/merge_bams.sh 82047
sbatch scripts/merge_bams.sh 82050

# Generate `gvcf` files
Here, we generate a separate `gvcf` per unique `sample * chromosome` combination in order to parallelize as much as possible.

First, generate a text file `prepared_bams.txt` listing the full path to all `bam` files:

```bash
ls OUTPUT/MERGED_BAMS/*.bam > prepared_bams.txt
ls OUTPUT/REMADE_BAMS/*.bam >> prepared_bams.txt
```

Then submit job array for [`build-gvcf.sh`](scripts/build-gvcf.sh). The 1st job generates `gvcf` for the first sample in `prepared_bams.txt`, and so  on, by using the `${SLURM_ARRAY_TASK_ID}` variable assigned to each job in the array.

We have 154 unique bam files, thus submit an array 1-154:
```bash
sbatch --array=1-154%10 scripts/build-gvcf.sh
sbatch --array=24 scripts/build-gvcf.sh

```

# Combine chromosomes per sample
Next, we combine all `gvcf` pertaining to a single sample using `gatk GatherVcfs`.

```bash
ls OUTPUT/gvcf/${sample}_*.g.vcf.gz
gatk GatherVcfs \
    -I ${sample}.list \
    -O ${sample}.g.vcf.gz \
    -R ${REF} \
    --CREATE_INDEX True

```

# Combine all samples
Next, we combine all samples using `GenomicsDBImport` (which is preferred to `CombineGVCFs`).
```bash
gatk GenomicsDBImport \
    -V data/gvcfs/mother.g.vcf \
    -V data/gvcfs/father.g.vcf \
    -V data/gvcfs/son.g.vcf \
    --genomicsdb-workspace-path my_database \
    --intervals chr20,chr21
```
# Joint genotyping
`GenotypeGVCFs `

# Extract alt contigs only
```bash
module load samtools
inbam='OUTPUT/REMADE_BAMS/81940.bam'
id=$(basename ${inbam%.bam})

```
```bash
module load GATK


cat alt_contigs.list | tr "\n" " " | xargs samtools view -hb ${inbam} | \
gatk AddOrReplaceReadGroups \
    I=/dev/stdin \
    O=${id}_alt_contigs_rgs.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${ID} \
    CREATE_INDEX=True 
```

```bash
module purge
module load samtools
module load GATK

REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
inbam='OUTPUT/REMADE_BAMS/81940.bam'
id=$(basename ${inbam%.bam})
contig='alt_contigs'

cat alt_contigs.list | tr "\n" " " | xargs samtools view -hb ${inbam} | \
gatk AddOrReplaceReadGroups \
    I=/dev/stdin \
    O=test.out.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${id} \
    CREATE_INDEX=true
    

gatk HaplotypeCaller \
    -R ${REF} \
    -L alt_contigs.list \
    -I test.out.bam \
    -O ${id}_${contig}.g.vcf.gz \
    -ERC GVCF
```

```bash
contig='chr1'

samtools view -hb ${inbam} ${contig} | \
gatk AddOrReplaceReadGroups \
    I=/dev/stdin \
    O=${id}.${contig}.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${id} \
    CREATE_INDEX=true
    

gatk HaplotypeCaller \
    -R ${REF} \
    -L ${contig} \
    -I ${id}.${contig}.bam \
    -O ${id}.${contig}.g.vcf.gz \
    -ERC GVCF
```

```bash
# Set up list of contigs for HaplotypeCaller
chrs=($(seq 22 -1 1))
chrs=${chrs[@]/#/chr}
intervals="chrX chrY chrM ${chrs}"
intervals=(${intervals})

for interval in ${intervals[@]}; do
    while read bamfile; do
        id=$(basename ${bamfile%.bam})
        id=${id%.merged}
        if [ ! -f "OUTPUT/gvcf/${id}_${interval}.g.vcf.gz.tbi" ]; then
            echo $bamfile $interval
        fi
        #echo $(realpath $bam) $interval
    done < prepared_bams.txt
done > for_haplotype_caller.txt

```



sbatch --array=1 scripts/make_examples.sh chrY

scripts/make_examples.sh chrY

```bash
jid='8307418'
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr8)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr9)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr10)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr11)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr12)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr13)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr14)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr15)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr16)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr17)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr18)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr19)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr20)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr21)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chr22)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chrX)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chrY)
jid=$(sbatch --array=1-154%10 --dependency afterany:${jid} scripts/make_examples.sh chrM)
```

```bash
while read chr len; do
    echo -e "$chr\t0\t$len" > hg38_$chr.bed
done < <(samtools idxstats OUTPUT/MERGED_BAMS/81988.merged.bam | \
    grep -v 'chrUn' | \
    grep -v 'random' | \
    grep -v 'EBV' | \
    grep -v '*' | cut -f 1,2)
```

```bash
sbatch --array=142 scripts/make_examples.sh chr7
8658188
```

sbatch scripts/glnexus_bcf.sh chr11
sbatch scripts/glnexus_bcf.sh chr12
sbatch scripts/glnexus_bcf.sh chr14
sbatch scripts/glnexus_bcf.sh chr15

grep -n $(basename $(grep -v -f <(ls OUTPUT/gvcf/chr13/ | cut -d '_' -f 1) prepared_bams.txt) | cut -d '.' -f 1) prepared_bams.txt
sbatch --array=129 scripts/make_examples.sh chr13
8658188_142
samtools view -bh <file.bam> <read>  >  file_read.bam 

module load samtools
samtools view  -bh /vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/OUTPUT/REMADE_BAMS/82060.bam chr6 > 82060_chr6.bam

samtools view  -f 256 /vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/OUTPUT/REMADE_BAMS/82060.bam > reads_256.bam

samtools view -bh in.bam chr1 > in_chr1.bam