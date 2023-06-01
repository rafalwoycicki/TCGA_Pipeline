nextflow.enable.dsl=2

process MapSequences {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    val index_prefix
    path fasta

    output:
    path("aligned.sam"), emit: aligned_sam

    script:
    """
    bowtie2 -L 20 -N 1 -x ../../../${index_prefix} -f ${fasta} -S aligned.sam
    """
}

process FilterSequneces1Mismatch {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path aligned_sam

    output:
    path("filtered.sam"), emit: filtered_sam

    script:
    """
    gawk '{if (\$0 ~ /^@/ || \$0 ~ /NM:i:[0-1][^0-9]/) print \$0}' ${aligned_sam} > filtered.sam
    """
}

process ExtractInfo {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path filtered_sam

    output:
    path("sorted.bam"), emit: sorted_bam

    script:
    """
    samtools view -S -b ${filtered_sam} > filtered.bam
    samtools sort filtered.bam -o sorted.bam
    """
}

process ExtractGFF {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path gffzipped

    output:
    path("reference.gff3"), emit: reference_gff3

    script:
    """
    gunzip -c ${gffzipped} > reference.gff3
    """
}

process GetGeneAnnotations {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path reference_gff3

    output:
    path("genes.bed"), emit: genes_bed

    script:
    """
    cat ${reference_gff3} | grep -v '#' | gawk '\$3 ~ "gene"' | gawk -v OFS="\t" '{split(\$9,a,";"); split(a[1],b,"ID=gene:"); split(a[2],c,"Name="); print \$1,\$4,\$5,\$7,b[2]";"c[2]}' > genes.bed
    """
}

process AnnotateGenes {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path sorted_bam
    path genes_bed

    output:
    path("annotated.bed"), emit: annotated_bed

    script:
    """
    bedtools bamtobed -i ${sorted_bam} > sorted.bed
    bedtools intersect -wa -wb -a sorted.bed -b ${genes_bed} > annotated.bed
    """
}

process CompareGeneNames {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path annotated_bed

    output:
    path("compared_genes.txt"), emit: compared_genes

    script:
    """
    gawk -v FS="\t" -v OFS="\t" '{split(\$4,a,"|"); split(\$11,b,";");print a[3], b[2], b[1]}' ${annotated_bed} > compared_genes.txt
    """
}

process ExtractGeneIDs {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path compared_genes

    output:
    path("gene_ids.txt"), emit: gene_ids_file

    script:
    """
    cut -f3 ${compared_genes} | sort | uniq > gene_ids.txt
    """
}

process RetrieveExpression {
    container 'tcga_pipeline'
    publishDir 'results', mode: 'copy'

    input:
    path expression_script
    val samples
    path gene_ids_file

    output:
    path("expression_matrix.txt"), emit: expression_matrix

    script:
    """
    Rscript ${expression_script} ${samples} ${gene_ids_file}
    """
}


workflow {
    fasta = file(params.fasta)
    gffzipped = file(params.gffzipped)
    expression_script = file(params.expression_script)
    
    MapSequences(params.index_prefix, fasta)
    FilterSequneces1Mismatch(MapSequences.out.aligned_sam)
    ExtractInfo(FilterSequneces1Mismatch.out.filtered_sam)
    ExtractGFF(gffzipped)
    GetGeneAnnotations(ExtractGFF.out.reference_gff3)
    AnnotateGenes(ExtractInfo.out.sorted_bam, GetGeneAnnotations.out.genes_bed)
    CompareGeneNames(AnnotateGenes.out.annotated_bed)
    ExtractGeneIDs(CompareGeneNames.out.compared_genes)
    RetrieveExpression(expression_script, params.samples, ExtractGeneIDs.out.gene_ids_file)
}

