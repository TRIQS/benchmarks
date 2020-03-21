#!/bin/sh
#BSUB -n 24
#BSUB -o test
#BSUB -J test

export DCORE_MPIRUN_COMMAND="/usr/share/lava/1.0/linux2.6-glibc2.12-x86_64/bin/intelmpi-mpirun -np 24"
export DCORE_ALPSCTHYB_TIMELIMIT=120

./alps_cthyb > output-alpscthyb
./pyed > output-pyed
