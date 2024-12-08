module load SpectrumMPI

make

mpisubmit.pl -p 1 program -- 1 1 1 128
mpisubmit.pl -p 2 program -- 1 1 1 128
mpisubmit.pl -p 4 program -- 1 1 1 128
mpisubmit.pl -p 8 program -- 1 1 1 128

или

launch.sh