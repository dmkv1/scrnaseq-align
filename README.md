# Pipeline for alignment of 5' scRNAseq samples

## Gene expression alingment and count - STARsolo

### Genome indexing

Inputs are defined in `nextflow.config`:

```conf
genomeFastaFiles = "path/to/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
sjdbGTFfile = "path/to//refdata-gex-GRCh38-2024-A/genes/genes.gtf"
sjdbOverhang = 149
```

Note that `genes.gtf` must be ungzipped before using it.

### Alingment process

Inputs are defined in `samples.csv`. Its columns:

```
sample_id,gex_fq1,gex_fq2
```

## VDJ alingment and count - cellranger

