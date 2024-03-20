#!/bin/sh                                                              
#$ -N GaAs_Ph
#$ -pe smp 20
#$ -q all.q
#$ -cwd                                                                                                                                                         
#$ -V                                                                                                                                                           
#$ -S /bin/bash                                                                                                                                                 
 
lumin -ph >& lumin.out

