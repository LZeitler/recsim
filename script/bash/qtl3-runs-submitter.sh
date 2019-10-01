#!/bin/bash

wdirb=$SCRATCH/slim/$(date +%H%M%S%d%m%y)

wdir=$wdirb.26
mkdir -p $wdir
cd $wdir
bsub -J "qtl326[1-100]" -n 1 -W 24:00 -R "rusage[mem=4000]" -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl3-runs-26.sh $wdir'

wdir=$wdirb.42
mkdir -p $wdir
cd $wdir
bsub -J "qtl342[1-100]" -n 1 -W 24:00 -R "rusage[mem=4000]" -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl3-runs-42.sh $wdir'

wdir=$wdirb.34
mkdir -p $wdir
cd $wdir
bsub -J "qtl334[1-100]" -n 1 -W 24:00 -R "rusage[mem=4000]" -oo $HOME/logs/%J_%I.stdout -eo $HOME/logs/%J_%I.stderr '$HOME/pro/recsim/script/bash/qtl3-runs-34.sh $wdir'
