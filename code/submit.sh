#!/bin/sh
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=199:99:99
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify GPU queue
#PBS -q gpu
#
# nodes: number of nodes
#   ppn: number of processes per node
#  gpus: number of gpus per node
#  GPUs are in 'exclusive' mode by default, but 'shared' keyword sets them to shared mode.
#PBS -l nodes=2:ppn=4:gpus=4:exclusive
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N Trust But Verify
#
# specify email for notifications
#PBS -M kyleabeauchamp@gmail.com
#
# mail settings (one or more characters)
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -m n
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o myoutput

# Change to working directory used for job submission
cd $PBS_O_WORKDIR

# start spark workers
echo $CUDA_VISIBLE_DEVICES
hostname
#python ${HOME}/src/kyleabeauchamp/TrustButVerify/scripts/run_grid.py $PBS_ARRAYID
#mpirun -rmk pbs progname
