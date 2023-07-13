import pandas
import glob
import os

threads = 2

FASTQS = pandas.read_csv('fastq_samples.tsv', sep='\t', dtype=str)
BAMS = pandas.read_csv('bam_samples.tsv', sep='\t', dtype=str)

NEW_BAM_IDS =  BAMS.NEW_ID.tolist()
NEW_FASTQ_IDS = FASTQS.NEW_ID.tolist()


# sample_dict = dict(zip(samples.HBCC_SAMPLE_NAME, samples.CARD_SAMPLE_NAME))
LANES = ['001', '002', '003', '004']
# dict going from samplename : allfiles

hg38 = '/data/CARDPB/resources/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa'
hg19 = ''



rule all:
    input:
        'OUTPUT/HBCC.vcf.gz'

rule create_symlinks:
    output:
        touch('.create_symlinks')
    script:
        "scripts/create_symlinks.py"

rule download_hg19:
    output:
        'human_g1k_v37.fasta.gz'
    shell:
        '''
        wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        '''

rule index_ref:
    input: "{reference}"
    output: "{reference}.ann"
    envmodules:
        'bwa/0.7.17'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        bwa index {input}
        """

rule fastq_to_bam:
    input:
        '.create_symlinks',
        ref=f"{hg38}",
        ref_index=f"{hg38}.ann",
        r1="INPUT/FASTQ_FILES_RENAMED/{sample}_L{lane}_R1_001.fastq.gz",
        r2="INPUT/FASTQ_FILES_RENAMED/{sample}_L{lane}_R2_001.fastq.gz"
    output:
        "BAMS_unmerged/{sample}_L{lane}.bam"
    envmodules:
        'bwa/0.7.17'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools view -b - | samtools sort -o {output}
        """

rule merge_bams:
    input:
        expand("BAMS_unmerged/{{sample}}_L{lane}.bam", lane=LANES)
    output:
        "BAMS_merged/{sample}.bam"
    envmodules: 'samtools/1.15'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """




rule call_merged_haplotypes:
    input: "BAMS_merged/{sample}.bam"
    output: "gVCF_from_FASTQ/{sample}.g.vcf"
    envmodules: 'samtools/1.15'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """




rule bam_to_fastq:
    input:
        '.create_symlinks',
        "INPUT/BAM_FILES_RENAMED/{sample}.final.bam"
    output:
        r1="REMADE_FASTQ/{sample}.R1.fastq.gz",
        r2="REMADE_FASTQ/{sample}.R2.fastq.gz"
    envmodules: 'samtools/1.15'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """




rule remap_bams:
    input: 
        ref='human_g1k_v37.fasta.gz.ann',
        r1="REMADE_FASTQ/{sample}.R1.fastq.gz",
        r2="REMADE_FASTQ/{sample}.R2.fastq.gz"
    output:
        "REMADE_BAMS/{sample}.bam"
    envmodules:
        'bwa/0.7.17'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """



rule call_remade_bam_haplotypes:
    input:
        "REMADE_BAMS/{sample}.bam"
    output:
        "gVCF_from_BAM/{sample}.g.vcf"
    envmodules: 'samtools/1.15'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """




rule combine_gvcf:
    input:
        expand("gVCF_from_FASTQ/{sample}.g.vcf", sample=NEW_FASTQ_IDS[0]),
        #expand("gVCF_from_FASTQ/{sample}.g.vcf", sample=NEW_FASTQ_IDS),
        #expand("gVCF_from_BAM/{sample}.g.vcf", sample=NEW_BAM_IDS)
    output:
        'OUTPUT/HBCC.vcf.gz'
    envmodules: 'samtools/1.15'
    resources:
        mem_mb=32000
    threads: 8
    shell:
        """
        """

