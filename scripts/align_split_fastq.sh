#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 50G
#SBATCH --gres lscratch:300
#SBATCH --time 12:00:00

set -v

id=${1}
part=${2}
REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')


r1="/vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/REMADE_FASTQ_FROM_BAM/${id}_split/${id}_R1_001.part_00${part}.fq.gz"
r2="/vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/REMADE_FASTQ_FROM_BAM/${id}_split/${id}_R2_001.part_00${part}.fq.gz"



# modules
module purge
module load GATK/4.4
module load bwa/0.7.17
module load samtools/1.17


# Set up working and output directories
OUTDIR="${PWD}/OUTPUT/UNMERGED_BAMS/"
FINALFILE="${OUTDIR}/${id}_part_00${part}.bam"


# Exit if finel output already exists
if [ -f "${FINALFILE}" ]; then
    echo "INFO: Final output ${FINALFILE} already exists."
    echo "Exiting..."
    exit 0
fi

mkdir -p ${OUTDIR}
TMP="/lscratch/${SLURM_JOB_ID}/"
mkdir -p ${TMP}
cd ${TMP}
echo "INFO: Working in ${TMP}"

# Map sample and output BAM
echo "INFO: Running bwa:"
bwa mem -t 8 ${REF} ${r1} ${r2} | samtools sort -@ 8 -T ${TMP} -o ${id}_part_00${part}.bam -



# Checkpoint to ensure BAM exists
if [ ! -f "${id}_part_00${part}.bam" ]; then
    echo "ERROR: the following bam file does not exist"
    echo "       ${id}_part_00${part}.bam"
    echo "Exiting..."
    exit 1
fi





# Move to permanent dir
echo "Moving bam to ${OUTDIR}"
mv ${id}_part_00${part}.bam ${OUTDIR}


# Done
echo "Done"
cd
exit 0
