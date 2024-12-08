gcc -o parallel main.c -lm -std=c99 -fopenmp

mpisubmit.pl -p 1 -t 2 parallel -- 1 1 1 256
mpisubmit.pl -p 1 -t 4 parallel -- 1 1 1 256