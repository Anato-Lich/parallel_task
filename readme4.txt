module load SpectrumMPI

make

mpisubmit.pl -p 1 -g 2 program -- 1 1 1 128
mpisubmit.pl -p 2 -g 1 program -- 1 1 1 128

или

launch.sh
