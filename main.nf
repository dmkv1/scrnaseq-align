#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARSOLO_INDEX } from './modules/alignment'

workflow {
    STARSOLO_INDEX(
        params.STAR.genomeFastaFiles,
        params.STAR.sjdbGTFfile
    )
}
