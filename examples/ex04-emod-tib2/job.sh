#!/bin/sh                                                              
#$ -N TiB2
#$ -pe smp 20
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash

lumin -emod >& lumin.out
