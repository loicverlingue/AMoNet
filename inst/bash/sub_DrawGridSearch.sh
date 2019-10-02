#!/bin/bash

#PBS -N Draw
#PBS -q batch
#PBS -l walltime=1:00:00
#PBS -l mem=5gb
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -m ae

cd /data/tmp/lverling/AMoNet/

/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript DrawGridSearch.R $CLASS
