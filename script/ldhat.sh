## phase data using beagle 5
module load gcc/4.8.2 gdc java/1.8.0_101
java -Xmx8g -jar ~/programs/beagle/beagle.16May19.351.jar gt=DM_HM_BP_BI_merged_fi_dip_SNO_sc2.recode.vcf out=SNO.phased

## get vcf to ldhat format
vcftools --gzvcf SNO.phased.vcf.gz --chr 'scaffold_2' --phased --ldhat --out 'SNO.phased' # output SNO.phased.ldhat.sites and *.locs

## get likelihood lookup tables
module load gcc/4.8.2 gdc ldhat/2.2a 
cd ~/pro/recsim/data/ldhat/
# scp lookup table from from ldhat repo, then downsample:
lkgen -lk lk_n100_t0.01 -nseq 10
interval -seq ../300/filter/SNO.phased.ldhat.sites -loc ../300/filter/SNO.phased.ldhat.locs -lk lk_n10_t0.01
stat # analyse rates, params: 100000 discarded, rates.txt
# does not work yet no output
# new try
cd $SCRATCH; interval -seq ~/pro/recsim/data/300/filter/SNO.phased.ldhat.sites -loc ~/pro/recsim/data/300/filter/SNO.phased.ldhat.locs -lk ~/pro/recsim/data/ldhat/lk_n10_t0.01 -bpen 30 -samp 3000 -its 1000000


## segment a file with vcftools
vcftools --gzvcf ../300/filter/SNO.phased.vcf.gz --from-bp 5915 --to-bp 817770 --chr 'scaffold_2' --phased --ldhat --out 'SNO.segm.phased'
1
