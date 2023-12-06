#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --cpus-per-task 12
#SBATCH --mem 100G
#SBATCH --gres lscratch:500
#SBATCH --time 5:00:00
#SBATCH --partition norm

module purge


cd ~/HBCC-ILM-WGS
wd=$(realpath ${PWD})

CONTIG=${1}

# set input dir
GVCFDIR=$(realpath "OUTPUT/gvcf/${CONTIG}")

# set output dir
OUTDIR="OUTPUT/vcf/"
mkdir -p ${OUTDIR}
OUTDIR=$(realpath ${OUTDIR})

# change to working dir
[[ -d /lscratch/${SLURM_JOB_ID} ]] && cd /lscratch/${SLURM_JOB_ID} || exit 1

# load modules
module load glnexus

# joint calling
glnexus_cli \
-m 100 \
-t 12 \
--config DeepVariant \
--bed "${wd}/hg38_${CONTIG}.bed" \
${GVCFDIR}/*.g.vcf.gz > ${CONTIG}.bcf

echo "done generating bcf"

module load bcftools

bcftools view "${CONTIG}.bcf" | bgzip -@ 12 -c > ${CONTIG}.vcf.gz

cp ${CONTIG}.vcf.gz ${OUTDIR}


