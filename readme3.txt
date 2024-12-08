module load SpectrumMPI

mpicxx task3.cpp -o parallel

mpisubmit.pl -p 1 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 2 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 4 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 8 -t 4 parallel -- 1 1 1 128

или

launch.sh