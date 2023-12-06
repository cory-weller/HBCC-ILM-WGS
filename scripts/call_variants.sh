#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 6
#SBATCH --mem 40G
#SBATCH --gres lscratch:50,gpu:k80:1
#SBATCH --time 3:00:00
#SBATCH --partition gpu

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
export N_SHARDS=$(nproc)
ID=$(basename ${BAM%.bam})
export ID=${ID%.merged}
export GVCF_TFRECORDS="gvcf.tfrecord@${N_SHARDS}.gz"
export EXAMPLES="make_examples.tfrecord@${N_SHARDS}.gz"

OUTDIR=deepvariant/${ID}_${CONTIG}
mkdir -p ${OUTDIR}
export OUTDIR=$(realpath ${OUTDIR})


[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID} || exit 1

rm -rf  logs
mkdir -p logs

cp ${OUTDIR}/* .

# CALL VARIANTS
MODEL="/opt/models/wgs/model.ckpt"
CALL_VARIANTS_OUTPUT="call_variants_output.tfrecord.gz"


call_variants \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --examples ${EXAMPLES} \
    --checkpoint "${MODEL}" || exit 1


# Append call_variants_output.tfrecord.gz
cp call_variants_output.tfrecord.gz ${OUTDIR} && \
cd $wd && \
sbatch scripts/postprocess_variants.sh ${CONTIG} ${N}
