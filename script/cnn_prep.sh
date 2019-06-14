# prepare popmask files for filtering
head -20 DM_HM_BP_BI_merged_fi_dip.vcf | grep -o 'SNO\w*' > /home/leo/pro/300_analyses/annotations/popmask/SNO.txt

# filter populations
vcftools --vcf /mnt/evo_euler/data/300/filter/DM_HM_BP_BI_merged_fi_dip.vcf --keep /home/leo/pro/recsim/300_analyses/annotations/popmask/SNO.txt --recode --out DM_HM_BP_BI_merged_fi_dip_SNO

# select only 1 chr for test runs (sc2) and remove non polymorphic sites, hwe and no missing data
vcftools --vcf /mnt/evo_euler/data/300/filter/DM_HM_BP_BI_merged_fi_dip_SNO.recode.vcf --maf 0.01 --hwe 0.01 --max-missing 1 --chr 'scaffold_2' --recode --out DM_HM_BP_BI_merged_fi_dip_SNO_sc2

# -> now it can happen that the number of SNPs in 2 popualtions is very different!


## get SNPs that are good in 2 populations so results are comparable
head -20 DM_HM_BP_BI_merged_fi_dip.vcf | grep -o 'SNO\w*' > /home/leo/pro/300_analyses/annotations/popmask/SNO_SZI.txt
head -20 DM_HM_BP_BI_merged_fi_dip.vcf | grep -o 'SZI\w*' >> /home/leo/pro/300_analyses/annotations/popmask/SNO_SZI.txt
vcftools --vcf /mnt/evo_euler/data/300/filter/DM_HM_BP_BI_merged_fi_dip.vcf --keep /home/leo/pro/recsim/300_analyses/annotations/popmask/SNO_SZI.txt --recode --out DM_HM_BP_BI_merged_fi_dip_SNO_SZI
# do common filtering
vcftools --vcf DM_HM_BP_BI_merged_fi_dip_SNO_SZI.recode.vcf --maf 0.01 --max-missing 1 --chr 'scaffold_2' --recode --out DM_HM_BP_BI_merged_fi_dip_SNO_SZI_sc2
# split in 2 
vcftools --vcf DM_HM_BP_BI_merged_fi_dip_SNO_SZI_sc2.recode.vcf --keep /home/leo/pro/recsim/300_analyses/annotations/popmask/SNO.txt --recode --out DM_HM_BP_BI_merged_fi_dip_SNO_sc2
vcftools --vcf DM_HM_BP_BI_merged_fi_dip_SNO_SZI_sc2.recode.vcf --keep /home/leo/pro/recsim/300_analyses/annotations/popmask/SZI.txt --recode --out DM_HM_BP_BI_merged_fi_dip_SZI_sc2

