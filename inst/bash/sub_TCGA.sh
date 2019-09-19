#!/bin/bash

#PBS -N TCGATrain
#PBS -q batch
#PBS -l walltime=30:00:00
#PBS -l mem=20gb
#PBS -l nodes=1:ppn=4
#PBS -j oe
#PBS -o /data/tmp/lverling
#PBS -m ae

cd /data/tmp/lverling

/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript RunTCGAopt_HPC.R $CLASS
