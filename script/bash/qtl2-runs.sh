#!/bin/bash

cd $1


for i in {1..1000}; do
    $HOME/programs/SLiM_build/slim $HOME/pro/recsim/eidos/qtl2.E
done

