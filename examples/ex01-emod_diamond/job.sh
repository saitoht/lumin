#!/bin/sh                                                              
#$ -N diamond
#$ -pe smp 20
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash

lumin -emod >& lumin.out
