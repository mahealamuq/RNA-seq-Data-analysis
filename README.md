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
```text
RNA-seq FASTQ Files
          │
          ▼
     FastQC
          │
          ▼
       fastp
          │
          ▼
      HISAT2
          │
          ▼
   Sorted BAM Files
          │
          ▼
   featureCounts
          │
          ▼
       DESeq2
          │
          ▼
 Differentially
 Expressed Genes
          │
          ▼
      BioMart
          │
          ▼
 Gene Symbols &
 Functional Annotation
          |
          ▼
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

Make the script executable:

```bash
chmod +x RNA-seq_Breast_cancer.sh
```

Run the script:

```bash
./RNA-seq_Breast_cancer.sh
```

---

## Explanation of Each Step

### 1. System Update

```bash
sudo apt update -y
sudo apt upgrade -y
```

This updates Ubuntu package information and upgrades installed packages.

---

### 2. Software Installation

The script installs required tools such as:

| Tool | Purpose |
|---|---|
| FastQC | Checks sequencing read quality |
| fastp | Trims low-quality reads |
| HISAT2 | Aligns RNA-seq reads to the genome |
| SAMtools | Processes SAM/BAM files |
| BEDTools | Genomic interval processing |
| Subread/featureCounts | Counts reads per gene |
| deepTools | Creates BigWig files |
| DESeq2 | Differential expression analysis |
| BioMart | Converts Ensembl IDs to gene symbols |

---

### 3. Miniconda and Conda Environment

The script installs Miniconda if it is not already installed.

Then it creates a Conda environment:

```bash
rnaseq_env
```

This environment contains bioinformatics packages used in the RNA-seq workflow.

---

### 4. Project Folder Structure

The script creates this project structure:

```text
breast_cancer_project/
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

### 5. Reference Genome Download

The script downloads the **hg38 HISAT2 genome index**.

This index is required for mapping RNA-seq reads to the human genome.

---

### 6. RNA-seq Data Download

The pipeline downloads paired-end RNA-seq data:

| Sample Type | Replicate | Accession |
|---|---|---|
| Normal breast tissue | Replicate 1 | ERR358485 |
| Normal breast tissue | Replicate 2 | ERR358486 |
| MCF7 breast cancer | Replicate 1 | ERR358488 |
| MCF7 breast cancer | Replicate 2 | ERR358487 |

---

### 7. Quality Control

FastQC checks raw sequencing quality.

It produces HTML reports showing:

- Per-base quality
- GC content
- Adapter contamination
- Sequence duplication
- Overrepresented sequences

---

### 8. Read Trimming

fastp removes poor-quality bases and sequencing artifacts.

Output files are stored in:

```text
trimmed/
```

---

### 9. RNA-seq Alignment

HISAT2 aligns trimmed reads to the human reference genome.

Output SAM files are created first, then converted into BAM files.

---

### 10. BAM Processing

SAMtools converts SAM files into sorted BAM files.

Sorted BAM files are required for:

- gene counting
- visualization
- statistics
- BigWig generation

---

### 11. Alignment Statistics

SAMtools flagstat reports mapping statistics for each sample.

Output files:

```text
results/normal_rep1_stats.txt
results/normal_rep2_stats.txt
results/MCF7_rep1_stats.txt
results/MCF7_rep2_stats.txt
```

---

### 12. BigWig File Generation

BigWig files are generated using deepTools.

These files can be loaded into IGV to visualize gene expression patterns across the genome.

---

### 13. Gene Annotation File

The script downloads:

```text
GENCODE v43 hg38 annotation
```

This GTF file tells featureCounts where genes are located in the genome.

---

### 14. Gene Quantification

featureCounts counts how many reads map to each gene.

Output:

```text
FeatureCounts/gene_counts.txt
```

This count matrix is used by DESeq2.

---

### 15. Differential Expression Analysis

DESeq2 compares:

```text
MCF7 vs normal
```

The output includes:

| Column | Meaning |
|---|---|
| baseMean | Average normalized expression |
| log2FoldChange | Expression difference |
| lfcSE | Standard error |
| stat | Test statistic |
| pvalue | Raw p-value |
| padj | Adjusted p-value |

Main output:

```text
results/DESeq2_results.csv
```

---

### 16. Top 20 Genes

The script creates:

```text
results/Top20_overexpressed.csv
results/Top20_underexpressed.csv
```

These files contain the most strongly upregulated and downregulated genes.

---

### 17. BioMart Gene Annotation

BioMart converts Ensembl gene IDs into readable gene names.

Example:

```text
ENSG00000141510 → TP53
```

Final annotated output:

```text
results/top20_overexpressed_genes_annotated.csv
results/top20_underexpressed_genes_annotated.csv
```

---

## Final Output Files

| Output File | Description |
|---|---|
| `fastqc/` | FastQC and MultiQC reports |
| `trimmed/` | Trimmed FASTQ files |
| `bam/*.bam` | Sorted alignment files |
| `bam/*.bai` | BAM index files |
| `bam/*.bw` | BigWig visualization files |
| `FeatureCounts/gene_counts.txt` | Gene count matrix |
| `results/DESeq2_results.csv` | Full differential expression results |
| `results/Top20_overexpressed.csv` | Top 20 upregulated genes in MCF7 |
| `results/Top20_underexpressed.csv` | Top 20 downregulated genes in MCF7 |
| `results/top20_overexpressed_genes_annotated.csv` | Annotated overexpressed genes |
| `results/top20_underexpressed_genes_annotated.csv` | Annotated underexpressed genes |

---

## Biological Interpretation

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

---

## Next Analysis Steps

After this pipeline finishes, the next steps are:

1. Open `Top20_overexpressed.csv`
2. Check gene symbols in the annotated file
3. Search important genes in GeneCards
4. Run disease enrichment analysis using Enrichr
5. Create plots such as:
   - volcano plot
   - heatmap
   - bar plot of top genes
6. Discuss cancer-related genes in the report

---

## Repository Structure

```text
breast-cancer-rnaseq-analysis/
│
├── RNA-seq_Breast_cancer.sh
├── README.md
│
└── breast_cancer_project/
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
