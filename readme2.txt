module load SpectrumMPI

mpicxx task2.cpp -o parallel

mpisubmit.pl -p 1 parallel -- 1 1 1 128
mpisubmit.pl -p 2 parallel -- 1 1 1 128
mpisubmit.pl -p 4 parallel -- 1 1 1 128
mpisubmit.pl -p 8 parallel -- 1 1 1 128

или

launch.sh

