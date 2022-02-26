#!/bin/bash
#$ -q  free*,pub64,class,som
#$ -N Parallel_Job
#$ -pe openmp 8-16
#$ -m beas
#$ -R y
#$ -ckpt blcr

module load R/3.5.1

cd /dfs6/pub/mcfarlar/

Rscript --vanilla Parallel_Base_Script_GP.R $1