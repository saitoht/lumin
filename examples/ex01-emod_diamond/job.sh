#!/bin/sh                                                              
#$ -N diamond
#$ -pe smp 20
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash
exe=$HOME/lumin

python $exe/lumin.py -emod >& lumin.out

