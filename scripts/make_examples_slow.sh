#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 3
#SBATCH --mem 20G
#SBATCH --gres lscratch:100
#SBATCH --time 4:00:00
#SBATCH --partition norm

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

# MAKE EXAMPLES
seq 0 $((N_SHARDS-1)) | \
  parallel --halt 2 --joblog 'logs/log' --res 'logs' \
    make_examples \
    --mode calling \
    --channels insert_size \
    --sample_name "${ID}" \
    --ref "${REF}" \
    --reads "${BAM}" \
    --regions "${CONTIG}" \
    --examples ${EXAMPLES} \
    --gvcf ${GVCF_TFRECORDS} \
    --task {} \
|| { cp -r logs ${OUTDIR}; exit 1; }

echo "File sizes:"
du -sh *

cp -r logs ${OUTDIR}
cp gvcf.tfrecord*  ${OUTDIR}
cp make_examples.tfrecord* ${OUTDIR} && \
cd ${wd} && \
sbatch scripts/call_variants_slow.sh ${CONTIG} ${N}
