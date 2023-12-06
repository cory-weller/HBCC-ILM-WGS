#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 100G
#SBATCH --time 24:00:00
#SBATCH --partition quick,norm
#SBATCH --gres lscratch:200

set -v

module purge
module load GATK
module load samtools
export TMPDIR="/lscratch/${SLURM_JOB_ID}/"
mkdir -p ${TMPDIR} 



# BAM: path to bam file
# prepared_bams.txt is list of ALL prepared bam files
# slurm job array variable extracts Nth line for job-specific BAM
export ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p prepared_bams.txt)
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

ls OUTPUT/gvcf/${sample}_*.g.vcf.gz > ${TMPDIR}/chrs.list

cd ${TMPDIR}


gatk GatherVcfs \
    -I chrs.list \
    -O ${sample}.g.vcf.gz \
    -R ${REF} \
    --CREATE_INDEX True


# Create bam with new read groups, on /lscratch
# This is done once per .bam, which is fed into HaplotypeCaller
# using the scatter-gather approach on all chromosomes in parallel
echo "Generating BAM with new read groups on /lscratch"

gatk AddOrReplaceReadGroups \
    I=${BAM} \
    O=${TMPDIR}/${ID}.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${ID} \
    CREATE_INDEX=True && echo 'done'



if [ ! -f "${ID}.bam.bai" ]; then
    samtools index ${ID}.bam
fi

make_gvcf() {
    local interval=${1}
    # drop '.list' for $contig, if file provided instead of specific contig
    local contig=$(basename ${interval%.list})
    local gvcf="${OUTDIR}/${ID}_${contig}.g.vcf.gz"

    if [ -f "${gvcf}.tbi" ]; then
        echo "INFO: ${gvcf}.tbi already exists!"
        echo "      This chromosome is already complete and will be skipped."
    else
        # remove gvcf if it exists, because lack of .tbi indicates
        # the file is incomplete or corrupted.
        [ -f "${gvcf}" ] && rm ${gvcf}

        echo "Processing ${ID}.bam into ${gvcf} for contig ${contig}"
        gatk HaplotypeCaller \
        --java-options "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms48G -Xmx48G -XX:ParallelGCThreads=4" \
        -R ${REF} \
        -I ${ID}.bam \
        -L ${interval} \
        -O ${gvcf} \
        -ERC GVCF
    fi
    echo "done with ${ID}:${contig}"
}
export -f make_gvcf


# to RUN ALL:
parallel -j 4 make_gvcf ::: ${intervals[@]}

# to run first 4: (x, y, M, alt_contigs):
# parallel -j 4 make_gvcf ::: ${intervals[@]:0:4}

