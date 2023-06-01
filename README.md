# TCGA_Pipeline
Dockerized pipeline to retrieve gene expression of sgRNA targets from TCGA db.

# Documentation for TCGA Analysis Pipeline

## Introduction
This pipeline is designed to process and analyze gene expression data from The Cancer Genome Atlas (TCGA). The pipeline is implemented using Nextflow, a data-driven computational workflow framework. It utilizes Docker containers for managing software dependencies, ensuring reproducibility and portability.

## Prerequisites

- Docker
- Nextflow
- R

## Dependencies
This pipeline depends on several bioinformatics tools, including:

- Bowtie2
- SAMtools
- Gawk
- Bedtools
- R with Bioconductor and packages (TCGAbiolinks and GenomicFeatures)

## Overview

1. **MapSequences**: This process maps sequences using bowtie2. The sequences are read from a FASTA file and mapped against a reference genome.

2. **FilterSequneces1Mismatch**: This process filters out sequences from the mapped sequences that have more than 1 mismatch.

3. **ExtractInfo**: This process converts the filtered SAM file into a BAM file and sorts it.

4. **ExtractGFF**: This process extracts GFF3 annotation from a gzipped file.

5. **GetGeneAnnotations**: This process extracts gene annotations from the GFF3 file and creates a BED file.

6. **AnnotateGenes**: This process annotates the sorted BAM file with the gene annotations from the BED file.

7. **CompareGeneNames**: This process compares gene names between different sources.

8. **ExtractGeneIDs**: This process extracts gene IDs from the compared genes.

9. **RetrieveExpression**: This process retrieves gene expression data for the extracted gene IDs.

## Input

The pipeline requires multiple input files:

- FASTA file: Contains the sequences to be mapped.
- Index Prefix: Prefix of the index files for the reference genome.
- GFF3.gz file: Contains the gene annotations in GFF3 format.
- Expression Script: An R script to retrieve expression data.
- Samples: A list of samples to be analyzed.

## Output

The pipeline produces several output files:

- `aligned.sam`: Contains the sequences aligned to the reference genome.
- `filtered.sam`: Contains the sequences with 1 or fewer mismatches.
- `sorted.bam`: Contains the sorted sequences in BAM format.
- `reference.gff3`: Contains the gene annotations in GFF3 format.
- `genes.bed`: Contains the gene annotations in BED format.
- `annotated.bed`: Contains the sequences annotated with gene information.
- `compared_genes.txt`: Contains the compared gene names.
- `gene_ids.txt`: Contains the extracted gene IDs.
- `expression_matrix.txt`: Contains the gene expression data.

## Usage

1. **Clone the Repository**:
   ```
   git clone <repository-url>
   ```
   Replace `<repository-url>` with the URL of your GitHub repository.

2. **Navigate to the Repository Directory**:
   ```
   cd TCGA_Pipeline
   ```
   
3. **Docker Setup

This pipeline uses Docker to manage these dependencies. To create a Docker image for this pipeline, a Dockerfile is provided in the repository. Here's how you can build and use the Docker image.

**Build Docker Image**: Navigate to the directory containing the Dockerfile and run the following command to build the Docker image:

```bash
docker build -f Dockerfile -t tcga_pipeline .
```

This command will create a Docker image named 'tcga_pipeline' that includes all the software dependencies needed for the pipeline.

**Run Docker Image**: You can test the Docker image by running the following command:

```bash
docker run -it tcga_pipeline
```

This command starts a new Docker container using the 'tcga_pipeline' image and opens an interactive terminal inside the container.
   
4. **Run the Pipeline**:
   ```
   nextflow run main.nf --fasta library.fa --index_prefix grch38prim --gffzipped Homo_sapiens.GRCh38.109.gff3.gz --expression_script retrieve_expression.R --samples TCGA-A7-A13D-01A-13R-A12P-07,TCGA-E9-A1RH-11A-34R-A169-07
   ```
   Adjust the parameters as necessary.

## Results

The pipeline will output several files, including an expression matrix file. This file contains gene expression values for the specified TCGA samples. The rows represent genes, and the columns represent samples.

## Conclusion

This documentation should provide a solid foundation to understand the TCGA analysis pipeline. If you have any further questions or run into issues, please refer to the code or consider reaching out for support.

