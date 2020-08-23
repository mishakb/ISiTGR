#!/bin/bash

# SRUN OPTIONS
# NOTE: #SBATCH is a command, not a comment

# ----------------------------------------------------------------------
# SOME PARAMETERS

# The number of tasks
#SBATCH -n 4
# ----------------------------------------------------------------------
# OPTIONAL PARAMETERS

# The number of cores per task (used for multi-threaded application, default=1,
# prefer a power of 6, no more than the max number of cores per node)
# Change this when using openMP in camb for CosmoMC.  Only advantageous when calculating
# CMB and MPK
##SBATCH -c 6
# The number of tasks per node (n * c)
# Can uncomment this to force x number of tasks per node.
##SBATCH --tasks-per-node=1

# The number of nodes
# SLURM issues a warning if you specify more nodes than tasks (N > n).
# It's better to let slurm calculate it for you.
# IMPORTANT: MPI between nodes will slow down the program because of the
# heavy I/O. It's recommended to more cores and less nodes for our jobs.
# Can change below to run MPI jobs across multiple nodes
#SBATCH -N 1

# Setting the name of the error-file to 'job.myjobid.err' and
# the output-file to 'job.myjobid.out'
# The file paths below are relative to the directory from which you submitted
# Change to your preferences. 
# We added two folders "errors" and "out_files" to save the following files
#SBATCH --error=errors/%J.err --output=../out_files/%J.out

#Print  detailed  event  logging to error file
#SBATCH -v

#Give your job a name, so you can more easily identify which job is which
#SBATCH -J Chains1
# ----------------------------------------------------------------------
# VERY OPTIONAL PARAMETERS

# Account name (project ID) to run under
##SBATCH -A <account>

# The maximum allowed run time (D-HH:MM:SS)
#SBATCH --time=4-00:00:00

# If this job needs 4GB of memory per mpi-task (=mpi ranks, =cores)
##SBATCH --mem-per-cpu=4000
# ----------------------------------------------------------------------
# Run the MPI application

#uncomment the line below if you want to run with openMP as advised above.
#change number of threads to match -c line above.
#export OMP_NUM_THREADS=6

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

# Changing to director where my application is since I submit from scripts
cd ../

#cd ../
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --mpi=pmi2 cosmomc test_ISiTGR.ini

echo "Program finished with exit code $? at: `date`"

# ---------------------------------------------------
# Reference: http://www.hpc2n.umu.se/batchsystem/examples_scripts,
#            http://www.hpc2n.umu.se/slurm-submit-file-design
#            https://computing.llnl.gov/tutorials/linux_clusters/man/srun.txt
# ---------------------------------------------------

