#!/bin/bash
set -e

# Increase soft limit for open files which STAR can temporarily require
ulimit -n 65536

THREADS=$1
THREADS=$((THREADS))

SAMPLE=$2
R1=$3
R2=$4

CBWHITELIST=$5
GENOME=$6
SJDBGTF=$7

CBLEN=$8
UMILEN=$9
# Convert to integers
CBLEN=$((CBLEN))
UMILEN=$((UMILEN))

SOLO_STRAND=${10}

STAR \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 128000000000 \
    --runThreadN $THREADS \
    --genomeDir $GENOME \
    --sjdbGTFfile $SJDBGTF \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --soloBarcodeMate 1 \
    --clip5pNbases 39 0 \
    --soloType CB_UMI_Simple \
    --clipAdapterType CellRanger4 \
    --soloCBwhitelist $CBWHITELIST \
    --soloCBstart 1 \
    --soloCBlen $CBLEN \
    --soloUMIstart $((CBLEN+1)) \
    --soloUMIlen $UMILEN \
    --soloStrand $SOLO_STRAND \
    --soloBarcodeReadLength 0 \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter None \
    --outFilterScoreMin 30 \
    --soloFeatures Gene GeneFull Velocyto \
    --soloMultiMappers EM \
    --outMultimapperOrder Random \
    --outFilterMultimapNmax 10 \
    --outSAMmultNmax 1 \
    --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
    --outFilterType BySJout \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outReadsUnmapped Fastx
