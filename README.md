# Pipeline for alignment of 5' scRNAseq samples

## Inputs

Inputs are defined in `samples.csv`.

```
sample_id,gex_fq1,gex_fq2,vdj_b_fq1,vdj_b_fq2,vdj_t_fq1,vdj_t_fq2
sample1,path/to/S1_GEX_R1.fastq.gz,path/to/S1_GEX_R2.fastq.gz,path/to/S1_BCR_R1.fastq.gz,path/to/S1_BCR_R2.fastq.gz,path/to/S1_TCR_R1.fastq.gz,path/to/S1_TCR_R2.fastq.gz
sample2,path/to/S2_GEX_R1.fastq.gz,path/to/S2_GEX_R2.fastq.gz,path/to/S2_BCR_R1.fastq.gz,path/to/S2_BCR_R2.fastq.gz,,
sample3,path/to/S3_GEX_R1.fastq.gz,path/to/S3_GEX_R2.fastq.gz,,path/to/S3_BCR_R2.fastq.gz,path/to/S1_TCR_R1.fastq.gz,path/to/S3_TCR_R2.fastq.gz
sample4,path/to/S4_GEX_R1.fastq.gz,path/to/S4_GEX_R2.fastq.gz,,,,
```

Absent TCR and BCR FASTQs would not affect the sample queue. However, gene expression FASTQs must be present for each sample.

## Gene expression alignment and count - STARsolo

### Genome indexing

Reference file inputs are defined in the STAR section of `params` in `nextflow.config`:

```conf
genomeFastaFiles = "path/to/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
sjdbGTFfile = "path/to/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
sjdbOverhang = 149
```

Note that `genes.gtf` must be ungzipped before using it.

Optionally, a STAR-pre-indexed genome can be used:

```conf
STAR {
    use_prebuilt_index = true
    genome_index = 'path/to/index_directory'
    ...
```

### Alignment process

Second part of the STAR section of `params` specifies barcode whitelist and 10X chemistry used in the experiment.

Parameters specific for Chromium Next GEM 5’ v2 protocol:

```conf
whitelist = 'path/to/737K-august-2016.txt'
cell_barcode_length = 16
umi_length = 10
soloStrand = 'Forward'
```

More STAR parameters are listed in the STAR invocation script `bin/STARsolo_align_gex.sh`:

```bash
STAR \
    --genomeLoad NoSharedMemory \
    --genomeDir $GENOME \
    --sjdbGTFfile $SJDBGTF \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --soloCBwhitelist $CBWHITELIST \
    --soloType CB_UMI_Simple \
    --soloStrand $SOLO_STRAND \
    --soloCBstart 1 \
    --soloCBlen $CBLEN \
    --soloUMIstart $((CBLEN+1)) \
    --soloUMIlen $UMILEN \
    --soloBarcodeMate 1 \
    --soloBarcodeReadLength 0 \
    --clip5pNbases 39 0 \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloMultiMappers EM \
    --outFilterMultimapNmax 10 \
    --outMultimapperOrder Random \
    --soloFeatures Gene GeneFull \
    --soloCellFilter None \
    --outFilterType BySJout \
    --outFilterScoreMin 30 \
    --outSAMmultNmax 1 \
    --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outReadsUnmapped Fastx
```



## VDJ alignment and count - cellranger

`cellranger vdj` call is very simple:

```groovy
script:
    """
    mv -v ${fq1} ${sample_id}_S1_L001_R1_001.fastq.gz
    mv -v ${fq2} ${sample_id}_S1_L001_R2_001.fastq.gz
    
    cellranger vdj --id=${sample_id} \
         --reference=${params.cellranger.cellranger_vdj_reference} \
         --fastqs=.
    """
```

## Acknowledgements

todo
