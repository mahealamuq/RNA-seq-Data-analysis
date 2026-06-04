# RNS-seq-Data-Analysis
RNA-seq analysis pipeline for identifying differentially expressed genes in MCF7 breast cancer cells compared to normal breast tissue.
Complete RNA-seq analysis pipeline using FastQC, HISAT2, SAMtools, featureCounts, and DESeq2 for identifying differentially expressed genes in MCF7 breast cancer cells.

This repository contains a complete RNA-seq bioinformatics workflow for analyzing differential gene expression between MCF7 breast cancer cells and normal breast tissue. The pipeline includes sequencing quality control, read alignment using HISAT2, BAM processing with SAMtools, gene quantification using featureCounts, and differential expression analysis using DESeq2.

The workflow is designed for paired-end RNA-seq experiments comparing cancer and normal breast tissue samples.

# Table of contents 

- [Introduction](#introduction)
- [Workflow Overview](#workflow-overview)
- [Software Requirements](#software-requirements)
- [Installation](#installation)
- [Data Download](#data-download)
- [Running the Pipeline](#running-the-pipeline)
- [Gene Quantification](#gene-quantification)
- [Differential Expression Analysis](#differential-expression-analysis)
- [Results](#results)
- [References](#references)

# RNA-seq Analysis Pipeline

## Introduction

This repository contains a complete RNA-seq bioinformatics workflow for identifying differentially expressed genes using next-generation sequencing (NGS) data.

RNA sequencing (RNA-seq) is widely used to study:
- Gene expression profiling
- Differential gene expression
- Transcriptome analysis
- Cancer-associated genes
- Biological pathways

This workflow demonstrates a standard RNA-seq analysis pipeline using:
- FastQC for quality control
- HISAT2 for genome alignment
- SAMtools for BAM processing
- featureCounts for gene quantification
- DESeq2 for differential expression analysis

The pipeline compares:
- MCF7 breast cancer samples
- Normal breast tissue samples

to identify significantly overexpressed and underexpressed genes.

---

## Workflow Overview

**Pipeline Workflow**

```text
FASTQ files
   ↓
FastQC quality control
   ↓
fastp read trimming
   ↓
HISAT2 alignment to hg38 genome
   ↓
SAMtools BAM processing
   ↓
BigWig file generation
   ↓
featureCounts gene quantification
   ↓
DESeq2 differential expression analysis
   ↓
Gene annotation
   ↓
Enrichr disease enrichment analysis
      ↓
 IGV for visually check overexpress
 and underexpress gene   
```
# Breast Cancer RNA-seq Analysis Pipeline

A complete RNA-seq analysis workflow for identifying differentially expressed genes between MCF7 breast cancer cells and normal breast tissue using HISAT2, featureCounts, DESeq2, and BioMart.

---

## Overview

This project implements an end-to-end RNA-seq bioinformatics pipeline for breast cancer transcriptome analysis.

The workflow performs:

- RNA-seq data download
- Quality assessment
- Read trimming
- Genome alignment
- Gene quantification
- Differential expression analysis
- Gene annotation
- Genome browser visualization

The objective is to identify genes that are significantly overexpressed or underexpressed in MCF7 breast cancer cells compared with normal breast tissue.

---

## Biological Background

Breast cancer is one of the most common cancers worldwide. Changes in gene expression contribute to cancer development, progression, metastasis, and treatment response.

RNA sequencing (RNA-seq) enables transcriptome-wide measurement of gene expression levels and allows identification of:

- Cancer-associated genes
- Potential biomarkers
- Therapeutic targets
- Dysregulated pathways

This project compares MCF7 breast cancer cells with normal breast tissue samples.

---

**Main Steps**

1. Download raw RNA-seq data
2. Perform sequencing quality control
3. Trim low-quality reads
4. Align reads to reference genome
5. Convert SAM to sorted BAM
6. Quantify reads per gene
7. Identify differentially expressed genes
8. Generate ranked gene list

---

## How to Run the Script
Clone The repository:

```bash
git clone https://github.com/USERNAME/breast-cancer-rnaseq-analysis.git
cd breast-cancer-rnaseq-analysis
```

Make the script executable:

```bash
chmod +x RNA-seq_Breast_cancer.sh
```

Run the script:

```bash
./RNA-seq_Breast_cancer.sh
```
Run the downstream R analysis separately if needed:

```bash
Rscript R_scripts/rnaseq_downstream_analysis.R
```

## Explanation of Each Step


## Step 1: System Update

The pipeline first updates Ubuntu packages:

```bash
sudo apt update -y
sudo apt upgrade -y
```

This ensures the operating system has the latest package information before installing bioinformatics tools.

---

## Step 2: Install Required Software

The script installs basic command-line and bioinformatics tools:

```bash
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
```

### Tool Purpose

| Tool | Purpose |
|---|---|
| FastQC | Checks raw sequencing read quality |
| fastp | Trims low-quality reads and adapters |
| HISAT2 | Aligns RNA-seq reads to the human genome |
| SAMtools | Converts, sorts, and indexes BAM files |
| featureCounts | Counts reads mapped to genes |
| deepTools | Generates BigWig files |
| DESeq2 | Performs differential expression analysis |
| Enrichr | Performs disease enrichment analysis |

---

## Step 3: Install and Load Miniconda

The pipeline checks whether Miniconda is installed.

If Miniconda is not found, it installs it automatically:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
```

Then Conda is loaded:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
```

---

## Step 4: Create Conda Environment

A Conda environment is created for the RNA-seq analysis:

```bash
conda create -n rnaseq_env python=3.10 -y
conda activate rnaseq_env
```

This keeps all required bioinformatics tools in a separate environment.

---

## Step 5: Install Conda Bioinformatics Packages

The pipeline installs additional packages from Bioconda and Conda-forge:

```bash
conda install -c bioconda -c conda-forge -y \
multiqc \
deeptools \
bioconductor-deseq2 \
bioconductor-biomart
```

These packages are used for report generation, genome visualization, and downstream R analysis.

---

## Step 6: Create Project Directories

The pipeline creates a structured project folder:

```bash
mkdir -p RNA-seq_analysis_project/{index,raw_data,fastqc,trimmed,bam,results,FeatureCounts,R_scripts,IGV}
```

### Directory Structure

```text
RNA-seq_analysis_project/
│
├── index/
├── raw_data/
├── fastqc/
├── trimmed/
├── bam/
├── results/
├── FeatureCounts/
├── R_scripts/
└── IGV/
```

---

## Step 7: Download Human Reference Genome Index

The pipeline downloads the prebuilt HISAT2 index for the human hg38 genome:

```bash
wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz

tar -xzf hg38_genome.tar.gz
```

The genome index is required for mapping RNA-seq reads to the human genome.

---

## Step 8: Download RNA-seq FASTQ Data

The pipeline downloads paired-end RNA-seq FASTQ files from ENA.

### Normal Breast Samples

| Sample | Accession |
|---|---|
| Normal replicate 1 | ERR358485 |
| Normal replicate 2 | ERR358486 |

### MCF7 Breast Cancer Samples

| Sample | Accession |
|---|---|
| MCF7 replicate 1 | ERR358488 |
| MCF7 replicate 2 | ERR358487 |

The files are renamed for clarity:

```text
Normal_rep1_1.fastq.gz
Normal_rep1_2.fastq.gz
Normal_rep2_1.fastq.gz
Normal_rep2_2.fastq.gz
MCF7_rep1_1.fastq.gz
MCF7_rep1_2.fastq.gz
MCF7_rep2_1.fastq.gz
MCF7_rep2_2.fastq.gz
```

---

## Step 9: Quality Control with FastQC

FastQC checks the quality of raw sequencing reads:

```bash
fastqc raw_data/*.fastq.gz -o fastqc
```

FastQC reports help identify:

- Low-quality reads
- Adapter contamination
- GC content bias
- Overrepresented sequences
- Sequence duplication levels

The output files are saved in:

```text
fastqc/
```

---

## Step 10: Read Trimming with fastp

fastp trims low-quality bases and removes sequencing artifacts.

Example:

```bash
fastp \
-i raw_data/Normal_rep1_1.fastq.gz \
-I raw_data/Normal_rep1_2.fastq.gz \
-o trimmed/Normal_rep1_trimmed_1.fastq.gz \
-O trimmed/Normal_rep1_trimmed_2.fastq.gz
```

Trimmed reads are saved in:

```text
trimmed/
```

---

## Step 11: RNA-seq Alignment with HISAT2

Trimmed reads are aligned to the hg38 human genome using HISAT2.

Example:

```bash
hisat2 \
-x index/hg38/genome \
-1 trimmed/Normal_rep1_trimmed_1.fastq.gz \
-2 trimmed/Normal_rep1_trimmed_2.fastq.gz \
-S bam/normal_rep1.sam
```

HISAT2 produces SAM alignment files.

---

## Step 12: Convert SAM to Sorted BAM

SAM files are converted into sorted BAM files using SAMtools:

```bash
samtools view -bS bam/normal_rep1.sam | samtools sort -o bam/normal_rep1_sorted.bam
```

Sorted BAM files are required for:

- Read counting
- Visualization
- Alignment statistics
- BigWig generation

---

## Step 13: Index BAM Files

Each BAM file is indexed:

```bash
samtools index bam/normal_rep1_sorted.bam
```

This creates `.bai` index files, which are needed for fast genome browser visualization.

---

## Step 14: Alignment Statistics

SAMtools generates alignment statistics:

```bash
samtools flagstat bam/normal_rep1_sorted.bam > results/normal_rep1_stats.txt
```

These files show:

- Total reads
- Mapped reads
- Properly paired reads
- Alignment percentage

---

## Step 15: Generate BigWig Files

BigWig files are created using deepTools:

```bash
bamCoverage -b bam/normal_rep1_sorted.bam -o bam/normal_rep1.bw
```

BigWig files can be loaded into IGV to visualize RNA-seq signal across the genome.

---

## Step 16: Download GENCODE GTF Annotation

The pipeline downloads the GENCODE human gene annotation file:

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

gunzip gencode.v43.annotation.gtf.gz
```

The GTF file contains gene locations and is required for gene counting.

---

## Step 17: Gene Quantification with featureCounts

featureCounts counts how many aligned reads map to each gene:

```bash
featureCounts \
-T 4 \
-p \
-a FeatureCounts/hg38_gencode.v43.annotation.gtf \
-o FeatureCounts/gene_counts.txt \
bam/normal_rep1_sorted.bam \
bam/normal_rep2_sorted.bam \
bam/MCF7_rep1_sorted.bam \
bam/MCF7_rep2_sorted.bam
```

The main output is:

```text
FeatureCounts/gene_counts.txt
```

This count matrix is used for DESeq2 analysis.

---

## Step 18: Differential Expression Analysis with DESeq2

The R script reads the gene count matrix:

```r
counts <- read.table(
    "FeatureCounts/gene_counts.txt",
    header = TRUE,
    row.names = 1,
    skip = 1
)
```

Annotation columns are removed:

```r
counts <- counts[,6:ncol(counts)]
```

Sample condition information is created:

```r
coldata <- data.frame(
    condition = factor(
        c("normal", "normal", "MCF7", "MCF7")
    ),
    row.names = colnames(counts)
)
```

DESeq2 compares:

```text
MCF7 vs normal
```

The contrast is:

```r
contrast = c("condition", "MCF7", "normal")
```

This means:

```text
Positive log2FoldChange = higher expression in MCF7
Negative log2FoldChange = higher expression in normal tissue
```

---

## Step 19: Identify Significant Genes

Significant genes are selected using:

```r
padj < 0.05
```

```r
sig_genes <- subset(
    as.data.frame(res_sorted),
    padj < 0.05 & !is.na(padj)
)
```

The adjusted p-value controls the false discovery rate.

---

## Step 20: Select Top Overexpressed and Underexpressed Genes

Top overexpressed genes in MCF7 are selected by sorting from highest positive log2 fold change:

```r
top_overexpressed <- head(
    sig_genes[order(-sig_genes$log2FoldChange), ],
    100
)
```

Top underexpressed genes in MCF7 are selected by sorting from most negative log2 fold change:

```r
top_underexpressed <- head(
    sig_genes[order(sig_genes$log2FoldChange), ],
    100
)
```

### Interpretation

| log2FoldChange | Meaning |
|---|---|
| Positive | Gene is upregulated in MCF7 |
| Negative | Gene is downregulated in MCF7 |
| 0 | No expression difference |

---

## Step 21: Gene Annotation

Gene IDs are converted into readable gene symbols.

Example:

```text
ENSG00000141510 → TP53
```

The annotated gene files help with biological interpretation and disease analysis.

---

## Step 22: Disease Enrichment Analysis with Enrichr

Gene symbols from the top overexpressed and underexpressed genes are submitted to Enrichr.

The analysis uses disease databases such as:

```r
databases <- c(
    "DisGeNET",
    "Jensen_DISEASES"
)
```

Enrichr identifies diseases associated with the submitted gene list.

---

## Step 23: Save Disease Enrichment Results

Disease enrichment results are saved as CSV files containing:

```text
Disease, Adjusted_P_Value, Associated_Genes
```

Example output files:

```text
Overexpressed_DisGeNET_Significant_Diseases.csv
Overexpressed_Jensen_Significant_Diseases.csv
Underexpressed_DisGeNET_Significant_Diseases.csv
Underexpressed_Jensen_Significant_Diseases.csv
```

---

## Final Output Files

| Output File | Description |
|---|---|
| `fastqc/` | FastQC quality reports |
| `trimmed/` | Trimmed FASTQ files |
| `bam/*.bam` | Sorted BAM alignment files |
| `bam/*.bai` | BAM index files |
| `bam/*.bw` | BigWig visualization files |
| `FeatureCounts/gene_counts.txt` | Gene count matrix |
| `results/DESeq2_results.csv` | Full differential expression result |
| `results/Significant_genes.csv` | Significant genes with padj < 0.05 |
| `results/Top100_overexpressed.csv` | Top 100 genes upregulated in MCF7 |
| `results/Top100_underexpressed.csv` | Top 100 genes downregulated in MCF7 |
| `results/*Disease*.csv` | Disease enrichment results |

---

## Biological Interpretation

The main goal of this RNA-seq pipeline is to identify genes that are differentially expressed between MCF7 breast cancer cells and normal breast tissue.
Positive `log2FoldChange` means the gene is more highly expressed in MCF7 breast cancer cells.

Negative `log2FoldChange` means the gene is more highly expressed in normal breast tissue.

For example:

```text
log2FoldChange = 2
```

means the gene is approximately 4 times higher in MCF7.

```text
log2FoldChange = -2
```

means the gene is approximately 4 times lower in MCF7.

Genes with positive log2 fold change are more highly expressed in MCF7 cells and may be involved in cancer-related processes such as:

- Cell proliferation
- Tumour growth
- Metastasis
- Hormone response
- Apoptosis regulation

Genes with negative log2 fold change are more highly expressed in normal breast tissue and may represent genes reduced or silenced in cancer cells.

Disease enrichment analysis helps identify whether these genes are associated with breast cancer or other cancer-related diseases.

---
## References

**RNA-seq Differential Expression**

Love MI et al. (2014).  
Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.  
Genome Biology.  
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8

**HISAT2**

Kim D et al. (2019).  
Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype.  
Nature Biotechnology.  
https://www.nature.com/articles/s41587-019-0201-4

**featureCounts**

Liao Y et al. (2014).  
featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features.  
Bioinformatics.  
https://academic.oup.com/bioinformatics/article/30/7/923/232889

**FastQC**
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**SAMtools**
http://www.htslib.org

---

## Author

Mahe Alam

---

## License

MIT License
