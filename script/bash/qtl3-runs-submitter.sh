#!/bin/bash

wdir=$SCRATCH/slim/$(date +%H%M%d%m%y)/
mkdir -p $wdir
cd $wdir
bsub -J "qtl3[1-100]" -n 1 -W 12:00 -R "rusage[mem=8000]" -oo logs/%J_%I.stdout -eo logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl3-runs.sh $wdir'
