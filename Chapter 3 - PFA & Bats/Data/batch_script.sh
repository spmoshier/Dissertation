#!/usr/bin/env bash

 

#SBATCH --time=95:00:00

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --mem=16gb

#SBATCH --job-name=pfa

 

echo ----

echo Job start at `date`

 

module load julia

perl /home/moshier.12/Desktop/pfa/pfa_batchrun/pfa_newsimtree/script.pl > testout

 

echo ----

echo Job ended at `date`