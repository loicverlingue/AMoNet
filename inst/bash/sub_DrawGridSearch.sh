#!/bin/bash

#PBS -N Draw
#PBS -q batch
#PBS -l walltime=1:00:00
#PBS -l mem=5gb
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /data/tmp/lverling
#PBS -m ae

cd /data/tmp/lverling

/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript /data/tmp/lverling/R/PossiLandscape.R $CLASS