# Breast-cancer-rnaseq-analysis
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
FASTQ Files
     ↓
FastQC Quality Control
     ↓
Read Trimming (fastp)
     ↓
HISAT2 Alignment
     ↓
SAMtools BAM Processing
     ↓
featureCounts Quantification
     ↓
DESeq2 Differential Expression Analysis
```

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

## Software Requirements

**Operating System**

- Ubuntu 20.04+ recommended

**Required Software**

| Software | Purpose |
|---|---|
| FastQC | Read quality control |
| MultiQC | Aggregate QC reports |
| fastp | Read trimming |
| HISAT2 | RNA-seq alignment |
| SAMtools | BAM processing |
| featureCounts | Gene quantification |
| R / DESeq2 | Differential expression analysis |

---

## Installation

**Update Ubuntu**

```bash
sudo apt update
```

**Install Required Tools**

```bash
sudo apt install -y \
fastqc \
hisat2 \
samtools \
subread \
wget \
gzip \
python3-pip
```

**Install MultiQC**

```bash
pip3 install multiqc
```

**Install fastp**

```bash
sudo apt install fastp
```

---

## Data Download

**Create Project Directory**

```bash
mkdir -p rnaseq_project/{raw_data,fastqc,trimmed,index,bam,counts,results,genome}
cd rnaseq_project
```

**Download Reference Genome**

```bash
cd genome

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

gunzip hg19.fa.gz
```

**Build HISAT2 Index**

```bash
cd ../index

hisat2-build ../genome/hg19.fa hg19_index
```

**Download RNA-seq Data**

Example:
- Normal sample = ERR358486
- MCF7 sample = ERR358487

```bash
cd ../raw_data

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358486/ERR358486_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358486/ERR358486_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358487/ERR358487_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/ERR358487/ERR358487_2.fastq.gz
```

---

## Running the Pipeline

**Step 1: Quality Control**

```bash
cd ../fastqc

fastqc ../raw_data/*.fastq.gz -o .
```

Generate combined report:

```bash
multiqc .
```

---

**Step 2: Read Trimming**

```bash
cd ../trimmed

fastp \
-i ../raw_data/ERR358486_1.fastq.gz \
-I ../raw_data/ERR358486_2.fastq.gz \
-o normal_R1.trimmed.fastq.gz \
-O normal_R2.trimmed.fastq.gz
```

```bash
fastp \
-i ../raw_data/ERR358487_1.fastq.gz \
-I ../raw_data/ERR358487_2.fastq.gz \
-o mcf7_R1.trimmed.fastq.gz \
-O mcf7_R2.trimmed.fastq.gz
```

---

**Step 3: Align Reads with HISAT2**

```bash
cd ../bam

hisat2 \
-x ../index/hg19_index \
-1 ../trimmed/normal_R1.trimmed.fastq.gz \
-2 ../trimmed/normal_R2.trimmed.fastq.gz \
-S normal.sam
```

```bash
hisat2 \
-x ../index/hg19_index \
-1 ../trimmed/mcf7_R1.trimmed.fastq.gz \
-2 ../trimmed/mcf7_R2.trimmed.fastq.gz \
-S mcf7.sam
```

---

**Step 4: Convert SAM to BAM**

```bash
samtools view -bS normal.sam | samtools sort -o normal.sorted.bam

samtools view -bS mcf7.sam | samtools sort -o mcf7.sorted.bam
```

---

**Step 5: Index BAM Files**

```bash
samtools index normal.sorted.bam

samtools index mcf7.sorted.bam
```

---

**Step 6: Alignment Statistics**

```bash
samtools flagstat normal.sorted.bam

samtools flagstat mcf7.sorted.bam
```

---

## Gene Quantification

**Run featureCounts**

```bash
cd ../counts

featureCounts \
-p \
-a hg19_gencode.v19.annotation.gtf \
-o gene_counts.txt \
../bam/normal.sorted.bam \
../bam/mcf7.sorted.bam
```

**Important Output Files**

| File | Description |
|---|---|
| gene_counts.txt | Gene count matrix |
| summary file | Mapping statistics |

---

## Differential Expression Analysis

**Run DESeq2 in R**

```r
library(DESeq2)

counts <- read.delim(
  "counts/gene_counts.txt",
  comment.char = "#",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(counts)

countData <- counts[, c(7,8)]

rownames(countData) <- counts$Geneid

countData <- as.matrix(countData)

mode(countData) <- "integer"

sampleInfo <- data.frame(
  row.names = colnames(countData),
  condition = factor(c("Normal", "MCF7"))
)

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = sampleInfo,
  design = ~ condition
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)

log2FC <- log2(
  (norm_counts[,2] + 1) /
  (norm_counts[,1] + 1)
)

results_table <- data.frame(
  Geneid = rownames(norm_counts),
  Normal = norm_counts[,1],
  MCF7 = norm_counts[,2],
  log2FoldChange = log2FC
)

results_table <- results_table[
  order(abs(results_table$log2FoldChange), decreasing = TRUE),
]

write.csv(
  results_table,
  "Manual_Log2FC_Results.csv",
  row.names = FALSE
)

head(results_table)
```

---

**Filter Significant Genes**

```r
sig <- subset(resOrdered,
              padj < 0.05 &
              abs(log2FoldChange) > 1)

write.csv(as.data.frame(sig),
          "Significant_Genes.csv")
```

---

## Results

Expected outputs include:

- Quality control reports
- Trimmed FASTQ files
- Sorted BAM alignment files
- Gene count matrix
- Differential expression results
- Significant gene list

**Example Biological Interpretation**

- Positive log2FoldChange:
  - higher expression in MCF7 cells

- Negative log2FoldChange:
  - lower expression in MCF7 cells

- Low adjusted p-value:
  - statistically significant expression difference

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

**Repository Structure**

```text
rnaseq-analysis-pipeline/
│
├── raw_data/
├── fastqc/
├── trimmed/
├── genome/
├── index/
├── bam/
├── counts/
├── results/
├── scripts/
│
├── README.md
└── rnaseq_pipeline.sh
```

**Author**

Mahe Alam

**License**

MIT License
