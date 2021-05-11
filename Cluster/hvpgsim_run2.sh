#!/bin/bash
#$ -cwd -S /bin/bash -l mem=1G,time=02:00:00 -N hvpg_sim2 -j y -t 1:500

currind=${SGE_TASK_ID}

export R_LIBS_USER=/ifs/scratch/msph/biostat/gm2978/R_libs
module load R/4.0.3
scl enable devtoolset-8 bash

cd /ifs/scratch/msph/biostat/gm2978/HVPG_TRT/

Rscript --vanilla cluster_hvpg_did2.R ${currind}
