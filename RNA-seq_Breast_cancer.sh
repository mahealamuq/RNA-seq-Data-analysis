#!/bin/bash

# =========================================================
# Breast Cancer RNA-seq Analysis Pipeline
# Author: Mahe Alam
# =========================================================

set -e

echo "================================================="
echo " Breast Cancer RNA-seq & ChIP-seq Pipeline "
echo "================================================="

# =========================================================
# UPDATE SYSTEM
# =========================================================

echo "Updating Ubuntu packages..."

sudo apt update -y
sudo apt upgrade -y
# =========================================================
# INSTALL BASIC TOOLS
# =========================================================

echo "Installing required Ubuntu packages..."

sudo apt install -y \
wget \
curl \
unzip \
gzip \
tar \
git \
default-jre \
python3-pip \
build-essential \
firefox \
fastqc \
samtools \
bedtools \
bowtie2 \
hisat2 \
subread \
fastp \
sra-toolkit

# =========================================================
# INSTALL MINICONDA
# =========================================================

echo "Checking Miniconda installation..."

if [ ! -d "$HOME/miniconda3" ]; then

    echo "Installing Miniconda..."

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

fi

# =========================================================
# LOAD CONDA
# =========================================================

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# =========================================================
# CREATE CONDA ENVIRONMENT
# =========================================================

ENV_NAME="rnaseq_env"

echo "Creating Conda environment..."

if ! conda env list | grep -q "$ENV_NAME"; then

    conda create -n $ENV_NAME python=3.10 -y

fi

conda activate $ENV_NAME

# =========================================================
# INSTALL BIOINFORMATICS TOOLS
# =========================================================

echo "Installing bioinformatics packages..."
conda install -c bioconda -c conda-forge -y \
multiqc \
deeptools \
bioconductor-deseq2 \
bioconductor-biomart

# =========================================================
# INSTALL MACS2
# =========================================================

echo "Installing MACS2..."

if ! conda env list | grep -q "macs2_env"; then

    conda create -n macs2_env python=3.9 -y

fi
conda activate macs2_env
conda install -c bioconda -c conda-forge macs2=2.2.7.1 -y
conda install -c conda-forge numpy scipy -y

# Test installation

macs2 --version

# Return to main environment

conda deactivate

echo "Checking software versions..."
python --version
fastqc --version
hisat2 --version
samtools --version
bedtools --version
bamCoverage --version

# =========================================================

echo "Creating project directories..."

mkdir -p RNA-seq_analysis_project/{index,raw_data,fastqc,trimmed,bam,results,FeatureCounts,R_scripts,IGV}

cd breast_cancer_project

# =========================================================
# DOWNLOAD HUMAN GENOME
# =========================================================

echo "Downloading hg38 genome index..."

cd index

wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz

tar -xzf hg38_genome.tar.gz

cd ..

# =========================================================
# DOWNLOAD RNA-seq DATA
# =========================================================

echo "Downloading RNA-seq FASTQ files..."
cd raw_data

# =========================================================
# NORMAL BREAST SAMPLES
# =========================================================

# Normal Replicate 1

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358485/ERR358485_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358485/ERR358485_2.fastq.gz

mv ERR358485_1.fastq.gz Normal_rep1_1.fastq.gz
mv ERR358485_2.fastq.gz Normal_rep1_2.fastq.gz

# Normal Replicate 2

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358486/ERR358486_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358486/ERR358486_2.fastq.gz

mv ERR358486_1.fastq.gz Normal_rep2_1.fastq.gz
mv ERR358486_2.fastq.gz Normal_rep2_2.fastq.gz

# =========================================================
# MCF7 BREAST CANCER SAMPLES
# =========================================================
# MCF7 Replicate 1

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358488/ERR358488_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358488/ERR358488_2.fastq.gz

mv ERR358488_1.fastq.gz MCF7_rep1_1.fastq.gz
mv ERR358488_2.fastq.gz MCF7_rep1_2.fastq.gz

# MCF7 Replicate 2

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358487/ERR358487_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358487/ERR358487_2.fastq.gz

mv ERR358487_1.fastq.gz MCF7_rep2_1.fastq.gz
mv ERR358487_2.fastq.gz MCF7_rep2_2.fastq.gz

cd ..

# =========================================================
# FASTQC QUALITY CONTROL
# =========================================================

echo "Running FastQC..."

fastqc raw_data/*.fastq.gz -o fastqc   
echo "Running fastp trimming..."

fastp \
-i raw_data/Normal_rep1_1.fastq.gz \
-I raw_data/Normal_rep1_2.fastq.gz \
-o trimmed/Normal_rep1_trimmed_1.fastq.gz \
-O trimmed/Normal_rep1_trimmed_2.fastq.gz

fastp \
-i raw_data/Normal_rep2_1.fastq.gz \
-I raw_data/Normal_rep2_2.fastq.gz \
-o trimmed/Normal_rep2_trimmed_1.fastq.gz \
-O trimmed/Normal_rep2_trimmed_2.fastq.gz

fastp \
-i raw_data/MCF7_rep1_1.fastq.gz \
-I raw_data/MCF7_rep1_2.fastq.gz \
-o trimmed/MCF7_rep1_trimmed_1.fastq.gz \
-O trimmed/MCF7_rep1_trimmed_2.fastq.gz

fastp \
-i raw_data/MCF7_rep2_1.fastq.gz \
-I raw_data/MCF7_rep2_2.fastq.gz \
-o trimmed/MCF7_rep2_trimmed_1.fastq.gz \
-O trimmed/MCF7_rep2_trimmed_2.fastq.gz
# =========================================================
# RNA-seq ALIGNMENT
# =========================================================

echo "Running HISAT2 alignment..."

# Normal Replicate 1

hisat2 \
-x index/hg38/genome \
-1 trimmed/Normal_rep1_trimmed_1.fastq.gz \
-2 trimmed/Normal_rep1_trimmed_2.fastq.gz \
-S bam/normal_rep1.sam

# Normal Replicate 2

hisat2 \
-x index/hg38/genome \
-1 trimmed/Normal_rep2_trimmed_1.fastq.gz \
-2 trimmed/Normal_rep2_trimmed_2.fastq.gz \
-S bam/normal_rep2.sam

# MCF7 Replicate 1

hisat2 \
-x index/hg38/genome \
-1 trimmed/MCF7_rep1_trimmed_1.fastq.gz \
-2 trimmed/MCF7_rep1_trimmed_2.fastq.gz \
-S bam/MCF7_rep1.sam

# MCF7 Replicate 2

hisat2 \
-x index/hg38/genome \
-1 trimmed/MCF7_rep2_trimmed_1.fastq.gz \
-2 trimmed/MCF7_rep2_trimmed_2.fastq.gz \
-S bam/MCF7_rep2.sam

# =========================================================
# SAM TO SORTED BAM
# =========================================================

echo "Converting SAM to BAM..."

samtools view -bS bam/normal_rep1.sam | samtools sort -o bam/normal_rep1_sorted.bam
samtools view -bS bam/normal_rep2.sam | samtools sort -o bam/normal_rep2_sorted.bam
samtools view -bS bam/MCF7_rep1.sam | samtools sort -o bam/MCF7_rep1_sorted.bam
samtools view -bS bam/MCF7_rep2.sam | samtools sort -o bam/MCF7_rep2_sorted.bam

# =========================================================
# INDEX BAM FILES
# =========================================================

echo "Indexing BAM files..."

samtools index bam/normal_rep1_sorted.bam
samtools index bam/normal_rep2_sorted.bam
samtools index bam/MCF7_rep1_sorted.bam
samtools index bam/MCF7_rep2_sorted.bam

# =========================================================
# ALIGNMENT STATISTICS
# =========================================================

echo "Generating alignment statistics..."

samtools flagstat bam/normal_rep1_sorted.bam > results/normal_rep1_stats.txt
samtools flagstat bam/normal_rep2_sorted.bam > results/normal_rep2_stats.txt
samtools flagstat bam/MCF7_rep1_sorted.bam > results/MCF7_rep1_stats.txt
samtools flagstat bam/MCF7_rep2_sorted.bam > results/MCF7_rep2_stats.txt

# =========================================================
# GENERATE BIGWIG FILES
# =========================================================

echo "Generating BigWig coverage files..."

bamCoverage -b bam/normal_rep1_sorted.bam -o bam/normal_rep1.bw
bamCoverage -b bam/normal_rep2_sorted.bam -o bam/normal_rep2.bw
bamCoverage -b bam/MCF7_rep1_sorted.bam -o bam/MCF7_rep1.bw
bamCoverage -b bam/MCF7_rep2_sorted.bam -o bam/MCF7_rep2.bw

# =========================================================
# DOWNLOAD GTF FILE
# =========================================================

echo "Downloading GTF annotation..."

#cd FeatureCounts

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

gunzip gencode.v43.annotation.gtf.gz

mv gencode.v43.annotation.gtf hg38_gencode.v43.annotation.gtf

cd ..

# =========================================================
# FEATURECOUNTS
# =========================================================
echo "Running featureCounts..."

featureCounts \
-T 4 \
-p \
-a FeatureCounts/hg38_gencode.v43.annotation.gtf \
-o counts/gene_counts.txt \
bam/normal_rep1_sorted.bam \
bam/normal_rep2_sorted.bam \
bam/MCF7_rep1_sorted.bam \
bam/MCF7_rep2_sorted.bam

# =========================================================
# DOWNLOAD IGV
# =========================================================

echo "Downloading IGV..."

cd IGV

wget https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip

unzip IGV_Linux_2.17.4_WithJava.zip

cd ..

# =========================================================
# PIPELINE COMPLETE
# =========================================================

echo "================================================="
echo " Pipeline Completed Successfully "
echo "================================================="
echo "Results Directory: breast_cancer_project/results"
echo "BigWig Files: breast_cancer_project/bam"
echo "Gene Counts: breast_cancer_project/FeatureCounts"
echo "================================================="


# ==========================================
# RNA-seq Downstream Analysis
# DESeq2 + BioMart Annotation + Enrichr Disease Analysis
# ==========================================

echo "Running BioMart annotation..."

Rscript R_scripts/biomart_annotation.R

