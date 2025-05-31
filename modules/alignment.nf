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
     --runThreadN ${params.STAR.index_threads} \
     --runMode genomeGenerate \
     --genomeDir "STAR_genome_indexed" \
     --genomeFastaFiles "${genomeFastaFile}" \
     --sjdbGTFfile "${sjdbGTFfile}" \
     --sjdbOverhang ${params.STAR.sjdbOverhang}
    """
}

process STARSOLO_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/GEX", mode: 'copy'
    
    input:
    tuple val(sample_id), path(gex_fq1), path(gex_fq2)
    path(whitelist)
    path(genome)
    path(sjdbGTFfile)
    
    output:
    path "Solo.out/", emit: counts
    path "Aligned.sortedByCoord.out.bam", emit: bam
    path "Log.final.out", emit: final_log

    script:
    def bamSortRAM = (params.STAR.sortbam_memory_gb * 1024L * 1024L * 1024L) as Long
    """
    STARsolo_align_gex.sh ${params.STAR.align_threads} ${bamSortRAM} \
        ${sample_id} ${gex_fq1} ${gex_fq2} \
        ${whitelist} ${genome} ${sjdbGTFfile} \
        ${params.STAR.cell_barcode_length} ${params.STAR.umi_length} ${params.STAR.soloStrand} 
    """
}

process CELLRANGER_VDJ {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/VDJ_B", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)
    path(cellranger_vdj_reference)

    output:
    path "${sample_id}/outs/", emit: vdj_b

    script:
    """
    cellranger vdj --id=${sample_id} \
         --reference=${params.cellranger.cellranger_vdj_reference} \
         --fastqs=.
    """
}