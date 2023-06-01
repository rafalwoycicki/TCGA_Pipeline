# Load necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Split the input string into a vector of sample names
samples <- unlist(strsplit(args[1], split=","))
gene_ids_file <- args[2]

# Load gene IDs
gene_ids <- readLines(gene_ids_file)

# Define query parameters
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  barcode = samples)

# Download data
GDCdownload(query)

# Prepare expression matrix
expression_matrix <- GDCprepare(query)

# Extract the expression data from the SummarizedExperiment
expression_data <- assay(expression_matrix)

# Remove the dot and everything on the right of the dot in the row names
rownames(expression_data) <- sub("\\..*", "", rownames(expression_data))

# Filter expression matrix for the genes in gene_ids
expression_data_filtered <- expression_data[rownames(expression_data) %in% gene_ids, ]

# Convert to data frame
expression_data_filtered_df <- as.data.frame(expression_data_filtered)

# Write filtered expression matrix to a file
write.table(expression_data_filtered_df, file = "expression_matrix.txt", sep = "\t", quote = FALSE)
