#!/bin/bash
#SBATCH -J CWB                  # Job name
#SBATCH -o CWB.o%j              # Name of stdout output file (%j expands to jobId)
#SBATCH -e CWB.o%j              # Name of stderr output file(%j expands to jobId)
#SBATCH -p normal               # Submit to the 'normal' or 'development' queue
#SBATCH -N 4                    # Total number of nodes
#SBATCH -n 272                  # Total number of mpi tasks requested
#SBATCH -t 48:00:00             # Run time (hh:mm:ss)
#SBATCH --mail-user=megha.j456@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH -A Cluster-Wild-Bootstr

# load R module
module load Rstats/3.5.1

# call R code from RMPISNOW
ibrun RMPISNOW < ./4_run_sim_study.R
