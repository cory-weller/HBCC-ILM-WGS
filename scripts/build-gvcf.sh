#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 40G
#SBATCH --time 6:00:00
#SBATCH --partition norm
#SBATCH --gres lscratch:60

set -v

module purge
module load GATK
module load samtools
export TMPDIR="/lscratch/${SLURM_JOB_ID}/"
mkdir -p ${TMPDIR} 

# Set up list of contigs for HaplotypeCaller
chrs=($(seq 22 -1 1))
chrs=${chrs[@]/#/chr}
export ALT_CONTIGS=$(cat alt_contigs.list)

intervals="chrX chrY chrM ${chrs}"
intervals=(${intervals})

# BAM: path to bam file
# prepared_bams.txt is list of ALL prepared bam files
# slurm job array variable extracts Nth line for job-specific BAM
BAM=$(sed -n ${SLURM_ARRAY_TASK_ID}p for_haplotype_caller.txt | cut -d ' ' -f 1)
export BAM=$(realpath $BAM)
export CONTIG=$(sed -n ${SLURM_ARRAY_TASK_ID}p for_haplotype_caller.txt | cut -d ' ' -f 2)
export ID=$(basename ${BAM%.bam})   # Get ID from bam file 
export ID="${ID%.merged}"           # if bam was merged, trim extension
export REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
if [ ! -f "${REF%.fa}.dict" ]; then
    gatk CreateSequenceDictionary \
        -R ${REF} \
        -O ${REF%.fa}.dict
fi

echo "Sample ID: ${ID}"

OUTDIR='OUTPUT/gvcf'
mkdir -p ${OUTDIR}
export OUTDIR=$(realpath ${OUTDIR})
export GVCF="${OUTDIR}/${ID}_${CONTIG}.g.vcf.gz"

if [ ! -f "${BAM}.bai" ]; then
    samtools index ${BAM}
fi

cp alt_contigs.list ${TMPDIR} && cd ${TMPDIR}




make_contig_gvcf() {
    if [ -f "${GVCF}.tbi" ]; then
        echo "INFO: ${GVCF}.tbi already exists!"
        echo "      This chromosome is already complete and will be skipped."
        return
    fi

    # remove gvcf if it exists, because lack of .tbi indicates the file is incomplete or corrupted.
    [ -f "${GVCF}" ] && rm ${GVCF}

    echo "Processing ${ID}.bam into ${GVCF} for contig ${CONTIG}"

    # Remake contig with correct read groups and create index
    samtools view -hb ${BAM} ${CONTIG} | \
    gatk AddOrReplaceReadGroups \
        I=/dev/stdin \
        O=${ID}_${CONTIG}.bam \
        SORT_ORDER=coordinate \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=${ID} \
        CREATE_INDEX=true
        
    # Call haplotype and export gvcf
    gatk --java-options "-XX:ParallelGCThreads=8" \
        HaplotypeCaller \
        -R ${REF} \
        -L ${CONTIG} \
        -I ${ID}_${CONTIG}.bam \
        -O ${ID}_${CONTIG}.g.vcf.gz \
        -ERC GVCF

    # Move gvcf to permanent disk
    rm ${ID}_${CONTIG}.bam
    mv ${ID}_${CONTIG}.g.vcf.gz ${OUTDIR}

    
    echo "done with ${ID}:${CONTIG}"
}
export -f make_contig_gvcf 


make_alt_contig_gvcf() {
    local gvcf="${OUTDIR}/${ID}_${CONTIG}.g.vcf.gz"

    if [ -f "${GVCF}.tbi" ]; then
        echo "INFO: ${GVCF}.tbi already exists!"
        echo "      This gvcf is already complete and will be skipped."
        return
    fi

    # remove gvcf if it exists, because lack of .tbi indicates the file is incomplete or corrupted.
    [ -f "${GVCF}" ] && rm ${GVCF}

    echo "Processing ${ID}.bam into ${GVCF} for ${CONTIG}"

    # Remake contig with correct read groups and create index
    samtools view -hb ${BAM} ${ALT_CONTIGS} | \
    gatk AddOrReplaceReadGroups \
        I=/dev/stdin \
        O=${ID}_alt_contigs.bam \
        SORT_ORDER=coordinate \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=${ID} \
        CREATE_INDEX=true
        
    # Call haplotype and export gvcf
    gatk --java-options "-XX:ParallelGCThreads=8" \
        HaplotypeCaller \
        -R ${REF} \
        -L alt_contigs.list \
        -I ${ID}_alt_contigs.bam \
        -O ${ID}_alt_contigs.g.vcf.gz \
        -ERC GVCF

    # Move gvcf to permanent disk
    rm ${ID}_alt_contigs.bam ${ID}_alt_contigs.bam.bai
    mv ${ID}_alt_contigs.g.vcf.gz ${OUTDIR}

    
    echo "done with ${ID}:alt_contigs"
}
export -f make_alt_contig_gvcf 

make_alt_contig_gvcf


make_contig_gvcf

exit 0