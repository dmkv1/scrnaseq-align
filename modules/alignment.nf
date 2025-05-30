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