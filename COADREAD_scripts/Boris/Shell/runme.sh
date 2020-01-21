#!/bin/bash
#SBATCH --nodes 1-1
#SBATCH --ntasks 1 --cpus-per-task=1
#SBATCH --qos bbdefault
#SBATCH --time 1600:0
#SBATCH --mail-type None

set -e

module purge; module load bluebear
./software/gdc-client download -t gdc-user-token.2019-10-28T17_48_48.851Z.txt -m gdc_manifest.2019-10-28.txt

