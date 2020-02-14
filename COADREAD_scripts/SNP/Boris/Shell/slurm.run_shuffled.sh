#!/bin/bash
#SBATCH -J shuffled                  # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH --mem 2000                 # in MB
#SBATCH -t 0-3:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o joblogs/job_%A_%a.log        # Standard output
#SBATCH -e joblogs/job_%A_%a.log        # Standard error
#SBATCH --account=noyvertb-cruk-bioinformatics

mkdir -p joblogs/


module purge; module load bluebear
module load R/v3.4.1

j=${SLURM_ARRAY_TASK_ID}
echo "job $j"

cluster_file="clusters_MSS.txt"
genotypes_file="matrix_genotypes.final.subset.txt"

outputdir="shuffled"

mkdir -p ${outputdir}

	echo "=================================================="

	echo "$(date): Started shuffling"
	cluster_file_shuffled=${outputdir}/shuffled_${j}.clusters
	cat <(head -n 1 ${cluster_file} ) <(paste <(tail -n +2 ${cluster_file}  | cut -f 1) <(tail -n +2 ${cluster_file} | cut -f 2 | shuf)) > ${cluster_file_shuffled}

	echo "$(date): Started association analysis"
	association_file="${cluster_file_shuffled%clusters}association"
	time Rscript --vanilla ProcessGenotypeAssociation.R ${cluster_file_shuffled} ${genotypes_file} > ${association_file}

        echo "$(date): Started annotation"
       	annotation_file="${cluster_file_shuffled%clusters}annotated"
	time bash annotate.association_file.sh ${association_file} > ${annotation_file}

	echo "=================================================="


echo "$(date) FINISHED"
