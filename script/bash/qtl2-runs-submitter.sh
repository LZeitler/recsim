#!/bin/bash

wdir=$SCRATCH/slim/$(date +%H%M%d%m%y)/
mkdir -p $wdir
bsub -J "qtl2[1-100]" -n 1 -W 12:00 -R "rusage[mem=8000]" -oo logs/%J_%I.stdout -eo logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl2-runs.sh $wdir'
