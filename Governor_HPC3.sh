#!/bin/bash
#SBATCH --job-name=Governor
#SBATCH -p free
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --mail-type=fail,invalid_depend
#SBATCH --mail-user=mcfarlar@hs.uci.edu

cd /dfs6/pub/mcfarlar/

for i in {1..20}
do
	sbatch  RunParallelScript_HPC3.sh $i
done
