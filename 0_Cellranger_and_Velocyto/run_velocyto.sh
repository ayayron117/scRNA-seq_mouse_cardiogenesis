#!/bin/bash

repeats="/home/amohammed/DATA/core_lab/scRNAseq_Ren/mm10_rmsk.gtf"
transcriptome="/home/amohammed/DATA/references/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf"
cellranger_output="/home/amohammed/DATA/core_lab/scRNAseq_Ren/output/Wildtype"

velocyto run10x -@ 20 -@ 6000 \
                -m $repeats \
                $cellranger_output \
                $transcriptome

cellranger_output="/home/amohammed/DATA/core_lab/scRNAseq_Ren/output/Knockout"

velocyto run10x -@ 20 -@ 6000 \
                -m $repeats \
                $cellranger_output \
                $transcriptome
