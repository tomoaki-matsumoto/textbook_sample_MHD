#! /bin/sh
ulimit -s 819200000
date
make clean
make
rm DATA/*
./main
gnuplot plot1d.gp

