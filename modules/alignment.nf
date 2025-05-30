process STARSOLO_INDEX {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(genomeFastaFile)
    path(sjdbGTFfile)
    
    output:
    path "STAR_genome_indexed/", emit: genome
    
    script:
    """
    STAR \
     --runThreadN ${params.STAR.threads} \
     --runMode genomeGenerate \
     --genomeDir "STAR_genome_indexed" \
     --genomeFastaFiles "${genomeFastaFile}" \
     --sjdbGTFfile "${sjdbGTFfile}" \
     --sjdbOverhang ${params.STAR.sjdbOverhang}
    """
}

process STARSOLO_ALIGN {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(gex_fq1), path(gex_fq2)
    path(whitelist)
    path(genome)
    path(sjdbGTFfile)
    
    script:
    """
    STARsolo_align_gex.sh ${params.STAR.threads} \
        ${sample_id} ${gex_fq1} ${gex_fq2} \
        ${whitelist} ${genome} ${sjdbGTFfile} \
        ${params.STAR.cell_barcode_length} ${params.STAR.umi_length} ${params.STAR.soloStrand} 
    """
}