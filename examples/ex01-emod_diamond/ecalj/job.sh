#!/bin/sh                                                              
#$ -N diamond_prim
#$ -pe smp 20
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash

lumin -emod >& lumin.out
