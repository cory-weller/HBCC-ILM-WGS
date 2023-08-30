#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 24G
#SBATCH --gres lscratch:75
#SBATCH --time 16:00:00

REF=$(realpath '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa')
DATADIR="${PWD}/INPUT/FASTQ_FILES_RENAMED/"
FASTQS=($(ls ${DATADIR}/*_R1_001.fastq.gz))
IDS=($(basename -a ${FASTQS[@]%_R1_001.fastq.gz}))

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
r2=${DATADIR}/${id}_R2_001.fastq.gz

echo "ID: ${id}"
echo "r1: ${r1}"
echo "r2: ${r2}"


# Checkpoint to ensure reads exist
if [  -f "${r1}" ] && [ ! -f "${r2}" ]; then
    echo "ERROR: the following read files do not exist"
    echo "       ${r1}"
    echo "       ${r2}"
    echo "Exiting..."
    exit 1
fi

# modules
module purge
module load bwa/0.7.17
module load samtools/1.17


# Set up working and output directories
OUTDIR="${PWD}/OUTPUT/UNMERGED_BAMS/"
mkdir -p ${OUTDIR}
TMP="/lscratch/${SLURM_JOB_ID}"
mkdir -p ${TMP}
cd ${TMP}
echo "INFO: Working in ${TMP}"


# Map sample and output BAM
echo "INFO: Running bwa:"
echo "      bwa mem -t 8 ${REF} ${r1} ${r2} | samtools sort -@ 8 -T ${TMP} -o ${id}.bam -" 
bwa mem -t 8 ${REF} ${r1} ${r2} | samtools sort -@ 8 -T ${TMP} -o ${id}.bam -


# Checkpoint to ensure BAM exists
if [ ! -f "${id}.bam" ]; then
    echo "ERROR: the following bam file does not exist"
    echo "       ${id}.bam"
    echo "Exiting..."
    exit 1
fi


# Move to permanent dir
echo "Moving ${id}.bam to ${OUTDIR}"
mv ${id}.bam ${OUTDIR}


# Done
echo "Done"
cd
exit 0
