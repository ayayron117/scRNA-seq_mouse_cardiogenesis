#!/bin/bash

ref='/home/amohammed/DATA/references/cellranger/mouse/refdata-gex-mm10-2020-A'
fastqdir='/home/amohammed/DATA/core_lab/scRNAseq_Ren/fastqz'
outdir='/home/amohammed/DATA/core_lab/scRNAseq_Ren/output'
      
while read name
do
      cd $outdir
      cellranger count --id=$name --transcriptome=$ref --fastqs=$fastqdir --sample=$name
done<sample_names
