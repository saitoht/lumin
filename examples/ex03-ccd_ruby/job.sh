#!/bin/sh                                                              
#$ -N Al2O3_Cr
#$ -pe smp 20
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash

lumin -ccd >& lumin.out
