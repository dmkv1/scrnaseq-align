#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARSOLO_INDEX } from './modules/alignment'
include { STARSOLO_ALIGN } from './modules/alignment'

workflow {
    STARSOLO_INDEX(
        params.STAR.genomeFastaFiles,
        params.STAR.sjdbGTFfile,
    )

    ch_samples_gex = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample
            def gex_fq1 = file(row.gex_fq1)
            def gex_fq2 = file(row.gex_fq1)

            return [sample_id, gex_fq1, gex_fq2]
        }

    STARSOLO_ALIGN(
        ch_samples_gex,
        params.STAR.whitelist,
        STARSOLO_INDEX.out.genome,
        params.STAR.sjdbGTFfile,
    )
}
