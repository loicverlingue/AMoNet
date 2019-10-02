#!/bin/bash
INIT=$(qsub -v CLASS="--Param MeanWinit SdWinit --DIR /data/tmp/lverling/AMoNet/ --NameProj LUNG --GENESman EGFR MTOR --Interval 10 --NewNet T --iteration 10 --SelectMECA HALLMARK --organ lusc|luad|hnsc --alpha 0 --lambda 0 --MinStepsForward 10 --MiniBatch 64 --MinConnect 4 --nblayers 4" -t 1-4 sub_TCGA.sh)
echo $INIT

DRAW=$(qsub -v CLASS="--NameProj LUNG --DIR /data/tmp/lverling/AMoNet/ --Validation T" -W depend=afteranyarray:$INIT sub_DrawGridSearch.sh)
echo $DRAW

FIRST=$(qsub -v CLASS="--Param lambda alpha --DIR /data/tmp/lverling/AMoNet/ --NameProj LUNG --NewNet F --iteration 10 --organ lusc|luad|hnsc " -W depend=afteranyarray:$INIT -t 1-10 sub_TCGA.sh)
echo $FIRST

DRAW=$(qsub -v CLASS="--NameProj LUNG --DIR /data/tmp/lverling/AMoNet/ --Validation T" -W depend=afteranyarray:$FIRST sub_DrawGridSearch.sh)
echo $DRAW

SECOND=$(qsub -v CLASS="--Param learningrate MinStepsForward --DIR /data/tmp/lverling/AMoNet/ --NameProj LUNG --NewNet F --iteration 20 --organ lusc|luad|hnsc " -W depend=afteranyarray:$FIRST -t 1-10 sub_TCGA.sh)
echo $FIRST
