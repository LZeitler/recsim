#!/bin/bash

dir=1900080919 # 1808080919  1855080919  1900080919
filename=qtl3-phenotypes.txt

cd /cluster/scratch/zeitlerl/slim/$dir
mkdir -p ../all/$dir
for file in $(find . -name $filename); do
    out=$(dirname $file | sed 's/\.//g' | sed 's,\/,,g')
    cp $file ../all/$dir/$(basename -s .txt $file).$out.txt
done
