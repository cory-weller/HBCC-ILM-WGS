#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --gres lscratch:50
#SBATCH --time 1:00:00
#SBATCH --partition quick,norm

module purge
module load deepvariant/1.5.0
module load parallel

cd ~/HBCC-ILM-WGS
wd=$(realpath $PWD)

export CONTIG=${1}
N=${2}; [ ! -z "${SLURM_ARRAY_TASK_ID}" ] && N=${SLURM_ARRAY_TASK_ID}

export REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
BAM=$(sed -n ${N}p prepared_bams.txt)
export BAM=$(realpath ${BAM})
export N_SHARDS=6
ID=$(basename ${BAM%.bam})
export ID=${ID%.merged}
export GVCF_TFRECORDS="gvcf.tfrecord@${N_SHARDS}.gz"
export EXAMPLES="make_examples.tfrecord@${N_SHARDS}.gz"

OUTDIR=deepvariant/${ID}_${CONTIG}
GVCF_OUTDIR="OUTPUT/gvcf/${CONTIG}"
mkdir -p ${OUTDIR}
mkdir -p "OUTPUT/gvcf/${CONTIG}"
export OUTDIR=$(realpath ${OUTDIR})
export GVCF_OUTDIR=$(realpath ${GVCF_OUTDIR})

[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID} || exit 1

rm -rf  logs
mkdir -p logs

cp ${OUTDIR}/* .

# CALL VARIANTS
MODEL="/opt/models/wgs/model.ckpt"
CALL_VARIANTS_OUTPUT="call_variants_output.tfrecord.gz"
OUTPUT_GVCF="${ID}_${CONTIG}.g.vcf.gz"
OUTPUT_VCF="${ID}_${CONTIG}.vcf.gz"

postprocess_variants \
    --ref "${REF}" \
    --infile "${CALL_VARIANTS_OUTPUT}" \
    --outfile "${OUTPUT_VCF}" \
    --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
    --gvcf_outfile "${OUTPUT_GVCF}" || exit 1

cp ${OUTPUT_GVCF} ${GVCF_OUTDIR} && rm -rf ${OUTDIR}
