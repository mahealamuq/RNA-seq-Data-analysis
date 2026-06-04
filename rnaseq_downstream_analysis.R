# ==========================================
# Create DESeq2 + BioMart + Enrichr Script
# ==========================================
# =========================================================
# RNA-seq Downstream Analysis
# DESeq2 + BioMart Annotation + Enrichr Disease Analysis
# =========================================================

# -----------------------------
# Load/install packages
# -----------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("enrichR")
library(DESeq2)
library(biomaRt)
library(enrichR)

# -----------------------------
# Read featureCounts file
# -----------------------------

counts <- read.table(
    "FeatureCounts/gene_counts.txt",
    header = TRUE,
    row.names = 1,
    skip = 1
)

# Remove annotation columns
counts <- counts[,6:ncol(counts)]

# -----------------------------
# Sample information
# -----------------------------

coldata <- data.frame(
    condition = factor(
        c(
            "normal",
            "normal",
            "MCF7",
            "MCF7"
        )
    ),
    row.names = colnames(counts)
)

# -----------------------------
# DESeq2 analysis
# -----------------------------

dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = ~ condition
)

dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c(
        "condition",
        "MCF7",
        "normal"
    )
)

res_sorted <- res[order(res$padj), ]

write.csv(
    as.data.frame(res_sorted),
    "results/DESeq2_results.csv"
)

# -----------------------------
# Significant genes
# -----------------------------

sig_genes <- subset(
    as.data.frame(res_sorted),
    padj < 0.05 & !is.na(padj)
)

sig_genes$Gene_id <- rownames(sig_genes)

sig_genes$Geneid_clean <- sub(
    "\\..*",
    "",
    sig_genes$Gene_id
)

write.csv(
    sig_genes,
    "results/Significant_genes.csv",
    row.names = FALSE
)

# -----------------------------
# Top 100 overexpressed and underexpressed genes
# -----------------------------

top_overexpressed <- head(
    sig_genes[order(-sig_genes$log2FoldChange), ],
    100
)

top_underexpressed <- head(
    sig_genes[order(sig_genes$log2FoldChange), ],
    100
)

write.csv(
    top_overexpressed,
    "results/Top100_overexpressed.csv",
    row.names = FALSE
)

write.csv(
    top_underexpressed,
    "results/Top100_underexpressed.csv",
    row.names = FALSE
)


# =========================================================
# BioMart Gene Annotation
# =========================================================

mart <- useMart(  "ensembl", dataset = "hsapiens_gene_ensembl")

annotate_genes <- function(gene_table, output_file) {

    gene_info <- getBM(
        attributes = c(
            "ensembl_gene_id",
            "hgnc_symbol",
            "external_gene_name",
            "description"
        ),
        filters = "ensembl_gene_id",
        values = gene_table$Geneid_clean,
        mart = mart
    )

    merged <- merge(
        gene_table,
        gene_info,
        by.x = "Geneid_clean",
        by.y = "ensembl_gene_id"
    )

    write.csv(
        merged,
        output_file,
        row.names = FALSE
    )

    return(merged)
}

top100_over_annotated <- annotate_genes(
    top_overexpressed,
    "results/Top100_overexpressed_annotated.csv"
)

top100_under_annotated <- annotate_genes(
    top_underexpressed,
    "results/Top100_underexpressed_annotated.csv"
)

# =========================================================
# Enrichr Disease Enrichment Analysis
# =========================================================

clean_symbols <- function(x) {
    x <- unique(x)
    x <- x[x != "" & !is.na(x)]
    return(x)
}

over_symbols <- clean_symbols(top100_over_annotated$hgnc_symbol)
under_symbols <- clean_symbols(top100_under_annotated$hgnc_symbol)

write.table(
    over_symbols,
    "results/Top100_overexpressed_gene_symbols.txt",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)

write.table(
    under_symbols,
    "results/Top100_underexpressed_gene_symbols.txt",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)

databases <- c(
    "DisGeNET",
    "Jensen_DISEASES"
)

disease_results_over <- enrichr(
    over_symbols,
    databases = databases
)

disease_results_under <- enrichr(
    under_symbols,
    databases = databases
)

# -----------------------------
# Save disease result table
# -----------------------------

save_disease_results <- function(enrichr_table, output_file) {

    if (is.null(enrichr_table) || nrow(enrichr_table) == 0) {

        empty_table <- data.frame(
            Disease = character(),
            Adjusted_P_Value = numeric(),
            Associated_Genes = character()
        )

        write.csv(
            empty_table,
            output_file,
            row.names = FALSE
        )

    } else {

        disease_table <- enrichr_table[, c(
            "Term",
            "Adjusted.P.value",
            "Genes"
        )]

        colnames(disease_table) <- c(
            "Disease",
            "Adjusted_P_Value",
            "Associated_Genes"
        )

        disease_table <- disease_table[
            order(disease_table$Adjusted_P_Value),
        ]

        write.csv(
            disease_table,
            output_file,
            row.names = FALSE
        )
    }
}


# -----------------------------
# Save significant disease results only
# -----------------------------

sig_over_DisGeNET <- subset(
    disease_results_over$DisGeNET,
    Adjusted.P.value < 0.05
)

sig_over_Jensen <- subset(
    disease_results_over$Jensen_DISEASES,
    Adjusted.P.value < 0.05
)

sig_under_DisGeNET <- subset(
    disease_results_under$DisGeNET,
    Adjusted.P.value < 0.05
)

sig_under_Jensen <- subset(
    disease_results_under$Jensen_DISEASES,
    Adjusted.P.value < 0.05
)

save_disease_results(
    sig_over_DisGeNET,
    "results/Overexpressed_DisGeNET_Significant_Diseases.csv"
)

save_disease_results(
    sig_over_Jensen,
    "results/Overexpressed_Jensen_Significant_Diseases.csv"
)

save_disease_results(
    sig_under_DisGeNET,
    "results/Underexpressed_DisGeNET_Significant_Diseases.csv"
)

save_disease_results(
    sig_under_Jensen,
    "results/Underexpressed_Jensen_Significant_Diseases.csv"
)

# -----------------------------
# Print summary
# -----------------------------

cat("====================================\n")
cat("RNA-seq downstream analysis complete\n")
cat("====================================\n")
cat("Main output files:\n")
cat("results/DESeq2_results.csv\n")
cat("results/Top100_overexpressed_annotated.csv\n")
cat("results/Top100_underexpressed_annotated.csv\n")
cat("results/Overexpressed_DisGeNET_Significant_Diseases.csv\n")
cat("results/Overexpressed_Jensen_Significant_Diseases.csv\n")
cat("results/Underexpressed_DisGeNET_Significant_Diseases.csv\n")
cat("results/Underexpressed_Jensen_Significant_Diseases.csv\n")














cat("====================================\n")

