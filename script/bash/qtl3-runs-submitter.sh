#!/bin/bash

wdir=$SCRATCH/slim/$(date +%H%M%d%m%y)/
mkdir -p $wdir
cd $wdir
bsub -J "qtl3-1[1-100]" -n 1 -W 24:00 -R "rusage[mem=4000]" -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl3-runs.sh $wdir'
