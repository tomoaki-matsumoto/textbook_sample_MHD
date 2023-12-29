#! /bin/sh
ulimit -s 819200000
date
make clean
make

mkdir -p DATA
rm DATA/*

# for Open MP
export OMP_STACKSIZE=512000
export OMP_NUM_THREADS=4
time ./main
echo $OMP_NUM_THREADS
# gnuplot plot1d.gp

