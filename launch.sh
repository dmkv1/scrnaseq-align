#/bin/bash
set -e

nextflow run main.nf \
  --samples samples.csv \
  --outdir results \
  -resume \
  -bg > run.log 2>&1
