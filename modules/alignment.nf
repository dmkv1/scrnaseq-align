process STARSOLO_INDEX {
    publishDir "reference", mode: 'copy'
    
    input:
    path(genomeFastaFile)
    path(sjdbGTFfile)
    
    output:
    path "STAR_genome_index/", emit: genome
    
    script:
    """
    STAR \
     --runThreadN ${task.cpus} \
     --runMode genomeGenerate \
     --genomeDir "STAR_genome_index" \
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
    def bamSortRAM = task.memory.toBytes() * 0.8
    """
    STARsolo_align_gex.sh ${task.cpus} ${bamSortRAM} \
        ${sample_id} ${gex_fq1} ${gex_fq2} \
        ${whitelist} ${genome} ${sjdbGTFfile} \
        ${params.STAR.cell_barcode_length} ${params.STAR.umi_length} ${params.STAR.soloStrand} 
    """
}

process CELLRANGER_VDJ_B {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/VDJ_BCR", mode: 'copy'

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

process CELLRANGER_VDJ_T {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/VDJ_TCR", mode: 'copy'

    input:
    tuple val(sample_id), path(fq1), path(fq2)
    path(cellranger_vdj_reference)

    output:
    path "${sample_id}/outs/", emit: vdj_t

    script:
    """
    cellranger vdj --id=${sample_id} \
         --reference=${params.cellranger.cellranger_vdj_reference} \
         --fastqs=.
    """
}