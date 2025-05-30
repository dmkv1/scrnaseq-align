#!/bin/bash
set -e

# Increase soft limit for open files which STAR can temporarily require
ulimit -n 65536

THREADS=$1
THREADS=$((THREADS))
BAM_SORT_RAM=$2

SAMPLE=$3
R1=$4
R2=$5

CBWHITELIST=$6
GENOME=$7
SJDBGTF=$8

CBLEN=$9
UMILEN=${10}
# Convert to integers
CBLEN=$((CBLEN))
UMILEN=$((UMILEN))

SOLO_STRAND=${11}

STAR \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM $BAM_SORT_RAM \
    --runThreadN $THREADS \
    --genomeDir $GENOME \
    --sjdbGTFfile $SJDBGTF \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --soloBarcodeMate 1 \
    --clip5pNbases 39 0 \
    --soloType CB_UMI_Simple \
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
