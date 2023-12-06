#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 20G
#SBATCH --time 24:00:00
#SBATCH --partition quick,norm

set -v

module purge
module load GATK

IDS=($@)
export REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')

if [ ! -f "${REF}.dict" ]; then
    gatk CreateSequenceDictionary \
        -R ${REF} \
        -O ${REF%.fa}.dict
fi


echo "IDs to process: ${IDS[@]}"


export OUTDIR="OUTPUT/BAM_RG"
mkdir -p ${OUTDIR}

replace_RGs() {

    local bam=${1}
    local id=$(basename ${bam%.bam})
    local id="${id%.merged}"
    local gvcf="${OUTDIR}/${id}.g.vcf.gz"

    echo "Adding read groups"

    gatk AddOrReplaceReadGroups \
    I=${bam} \
    O=${OUTDIR}/${id}.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${id} \
    CREATE_INDEX=True
    echo "done with ${id}"
}

export -f replace_RGs


parallel -j 1 replace_RGs ::: ${IDS[@]}
