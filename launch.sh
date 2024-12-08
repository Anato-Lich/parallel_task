#! /usr/bin/bash
mpisubmit.pl -p 1 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 2 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 4 -t 4 parallel -- 1 1 1 128
mpisubmit.pl -p 8 -t 4 parallel -- 1 1 1 128

mpisubmit.pl -p 1 -t 4 parallel -- 1 1 1 256
mpisubmit.pl -p 2 -t 4 parallel -- 1 1 1 256
mpisubmit.pl -p 4 -t 4 parallel -- 1 1 1 256
mpisubmit.pl -p 8 -t 4 parallel -- 1 1 1 256

mpisubmit.pl -p 1 -t 4 parallel -- 1 1 1 512
mpisubmit.pl -p 2 -t 4 parallel -- 1 1 1 512
mpisubmit.pl -p 4 -t 4 parallel -- 1 1 1 512
mpisubmit.pl -p 8 -t 4 parallel -- 1 1 1 512

mpisubmit.pl -p 1 -t 4 parallel -- 3.14 3.14 3.14 128
mpisubmit.pl -p 2 -t 4 parallel -- 3.14 3.14 3.14 128
mpisubmit.pl -p 4 -t 4 parallel -- 3.14 3.14 3.14 128
mpisubmit.pl -p 8 -t 4 parallel -- 3.14 3.14 3.14 128

mpisubmit.pl -p 1 -t 4 parallel -- 3.14 3.14 3.14 256
mpisubmit.pl -p 2 -t 4 parallel -- 3.14 3.14 3.14 256
mpisubmit.pl -p 4 -t 4 parallel -- 3.14 3.14 3.14 256
mpisubmit.pl -p 8 -t 4 parallel -- 3.14 3.14 3.14 256

mpisubmit.pl -p 1 -t 4 parallel -- 3.14 3.14 3.14 512
mpisubmit.pl -p 2 -t 4 parallel -- 3.14 3.14 3.14 512
mpisubmit.pl -p 4 -t 4 parallel -- 3.14 3.14 3.14 512
mpisubmit.pl -p 8 -t 4 parallel -- 3.14 3.14 3.14 512