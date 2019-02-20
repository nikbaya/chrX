#!/bin/bash
#$ -l h_rt=16:10:30
#$ -j y
#$ -l h_vmem=32g
#$ -cwd
#$ -o /broad/finucanelab/ktashman/inrich_analyses/simulations/download.log
#$ -t 0-400something

source /broad/software/scripts/useuse
reuse -q .anaconda-5.0.1
reuse -q .python-2.7.14-sqlite3-rtrees
reuse -q .google-cloud-sdk
mapfile -t myArray < phenotypes.txt
this_url=${myArray[${SGE_TASK_ID}]}
gsutil -m cp gs://ukbb-gwas-imputed-v3-results/export1/${this_url}.suffix /path/on/cluster