#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STARSOLO_INDEX } from './modules/alignment'
include { STARSOLO_ALIGN } from './modules/alignment'
include { CELLRANGER_VDJ_B } from './modules/alignment'
include { CELLRANGER_VDJ_T } from './modules/alignment'

workflow {
    if (params.STAR.use_prebuilt_index) {
        ch_genome = Channel.fromPath(params.STAR.genome_index, type: 'dir', checkIfExists: true)
    }
    else {
        STARSOLO_INDEX(
            params.STAR.genomeFastaFiles,
            params.STAR.sjdbGTFfile,
        )
        ch_genome = STARSOLO_INDEX.out.genome
    }

    ch_samples_gex = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def gex_fq1 = file(row.gex_fq1)
            def gex_fq2 = file(row.gex_fq2)

            return [sample_id, gex_fq1, gex_fq2]
        }

    STARSOLO_ALIGN(
        ch_samples_gex,
        params.STAR.whitelist,
        ch_genome,
        params.STAR.sjdbGTFfile,
    )

    ch_samples_vdj_b = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def fq1 = row.vdj_b_fq1
            def fq2 = row.vdj_b_fq2

            return [sample_id, fq1, fq2]
        }
        .filter { sample_id, fq1, fq2 -> 
            fq1 && fq2 && fq1.trim() != '' && fq2.trim() != ''
        }
        .map { sample_id, fq1, fq2 ->
            [sample_id, file(fq1), file(fq2)]
        }

    CELLRANGER_VDJ_B(
        ch_samples_vdj_b,
        params.cellranger.cellranger_vdj_reference
    )

    ch_samples_vdj_t = Channel.fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def fq1 = row.vdj_t_fq1
            def fq2 = row.vdj_t_fq2

            return [sample_id, fq1, fq2]
        }
        .filter { sample_id, fq1, fq2 -> 
            fq1 && fq2 && fq1.trim() != '' && fq2.trim() != ''
        }
        .map { sample_id, fq1, fq2 ->
            [sample_id, file(fq1), file(fq2)]
        }

    CELLRANGER_VDJ_T(
        ch_samples_vdj_t,
        params.cellranger.cellranger_vdj_reference
    )
}
