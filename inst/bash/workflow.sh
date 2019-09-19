#!/bin/bash
INIT=$(qsub -v CLASS="--Param lambda learningrate --DIR /data/tmp/lverling/PossiLandscape/ --NameProj LUNG --GENESman EGFR MTOR --Interval 10 --NewNet T --iteration 40 --SelectMECA HALLMARK --organ lusc|luad|hnsc --alpha 2 --MinStepsForward 10 --MiniBatch 64 " -t 1-10 sub_TCGA.sh)
echo $INIT

DRAW=$(qsub -v CLASS="--NameProj LUNG --DIR /data/tmp/lverling/PossiLandscape/ --Validation T" -W depend=afteranyarray:$INIT sub_DrawGridSearch.sh)
echo $DRAW
