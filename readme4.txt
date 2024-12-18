module unload SpectrumMPI
module load openmpi
module load pgi

pgcc -acc -ta=tesla:cc60 -O3 -Minfo=accel -o gpu_task MPI_ACC.c -L /opt/open_mpi/lib -I /opt/open_mpi/include -l mpi -cudalibs

bsub < acc_job.lsf
