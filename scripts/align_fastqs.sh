#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --gres lscratch:200
#SBATCH --time 8:00:00

REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
DATADIR="${PWD}/INPUT/FASTQ_FILES_RENAMED/"
FASTQS=($(ls ${DATADIR}/*_R1_001.fastq.gz))
IDS=($(basename -a ${FASTQS[@]%_R1_001.fastq.gz}))

echo ${IDS[0]}
# Get job array ID (or manually provided while testing)
if [ ! -z "${SLURM_ARRAY_TASK_ID}" ]; then
    N=${SLURM_ARRAY_TASK_ID}
elif [ ! -z "${1}" ]; then
    N=${1}
else
    echo 'ERROR: No SLURM job array number or manually defined number!'
    exit 1
fi


# Extract job-specific sample
id=${IDS[${N}]}
r1=${DATADIR}/${id}_R1_001.fastq.gz
r2=${DATADIR}/${id}_R1_001.fastq.gz



# modules
module purge
module load bwa/0.7.17
module load samtools/1.17


# Set up working and output directories
OUTDIR="${PWD}/OUTPUT/UNMERGED_BAMS/"
mkdir -p ${OUTDIR}
TMPDIR="/lscratch/${SLURM_JOB_ID}"
mkdir -p ${TMPDIR}
cd ${TMPDIR}
echo "INFO: Working in ${TMPDIR}"



# Map sample and output BAM
echo "INFO: Running bwa:"
echo "      bwa mem -t 8 -o ${id}.sam ${REF} ${r1} ${r2}" 
bwa mem -t 8 -o ${id}.sam ${REF} ${r1} ${r2}


# Convert SAM to BAM
echo "INFO: converting to bam:"
echo "      samtools view -b ${id}.sam | samtools sort -o ${id}.bam -"
samtools view -b ${id}.sam | samtools sort -o ${id}.bam -


# Remove sam
rm ${id}.sam


# Move to permanent dir
echo "Moving ${id}.bam to ${OUTDIR}"
mv ${id}.bam ${OUTDIR}


# Done
echo "Done"
cd
exit 0
