#!/bin/bash

# dir=1900080919 # 1808080919  1855080919  1900080919  150cM runs
# dir=1516250919 #1520250919 1518250919 1516250919 75cM runs
# dir=191334250919 #1909250919 191316250919 191334250919 75cM runs new means
dir=185759011019.26 # 185759011019.42 185759011019.34
filename=qtl3-phenotypes.txt

cd /cluster/scratch/zeitlerl/slim/$dir
mkdir -p ../all/$dir
for file in $(find . -name $filename); do
    out=$(dirname $file | sed 's/\.//g' | sed 's,\/,,g')
    cp $file ../all/$dir/$(basename -s .txt $file).$out.txt
done
