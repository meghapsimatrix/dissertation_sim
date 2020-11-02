#!/bin/bash
#SBATCH -J CWB	    	 	# Job name
#SBATCH -o CWB.o%j 		  	# Name of stdout output file (%j expands to jobId)
#SBATCH -e CWB.o%j 		  	# Name of stderr output file(%j expands to jobId)
#SBATCH -p development		# Submit to the 'normal' or 'development' queue
#SBATCH -N 4				# Total number of nodes
#SBATCH -n 176 		   		# Total number of mpi tasks requested
#SBATCH -t 00:30:00  		# Run time (hh:mm:ss)
#SBATCH --mail-user=megha.j456@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# load R module
module load Rstats        

# call R code from RMPISNOW
ibrun RMPISNOW < ./cwb/4_run_sim_study.R
