


cd /media/elgon/Jesse/HerbariumSeq/

#make a version with few missing SNPs to improve inference for the low coverage accessions

#1. #########################

~/Desktop/plink --tfile  /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --geno 0.05 --hwe 1e-50 'midp' --allow-no-sex --mind 0.95 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_1001G.txt --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem



cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem.imiss ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/


##LD thinning
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem --indep-pairwise 50 10 0.1
#(output indicating number of SNPs targeted for inclusion/exclusion)
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem --extract /media/elgon/Jesse/HerbariumSeq/plink.prune.in --make-bed --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem_pruned


#genetic distance thinned
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_1001Gbadrem_pruned --allow-extra-chr --distance square 1-ibs flat-missing 


cp /media/elgon/Jesse/HerbariumSeq/plink.mdist ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_1001Gbadrem.mdist
cp /media/elgon/Jesse/HerbariumSeq/plink.mdist.id ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_1001Gbadrem.mdist.id


#2. ##############next step is to redo SNP filtering for all selected genotypes###########
cd /media/elgon/Jesse/HerbariumSeq/

~/Desktop/plink --tfile  /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --geno 0.05 --hwe 1e-50 'midp' --allow-no-sex --mind 0.95 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange.imiss ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/


##LD thinning
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange --indep-pairwise 50 10 0.1
#(output indicating number of SNPs targeted for inclusion/exclusion)
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange --extract /media/elgon/Jesse/HerbariumSeq/plink.prune.in --make-bed --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange_pruned


#3. genetic distance, MDS, Admixture for all genotypes (global)

#genetic distance thinned
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange_pruned --allow-extra-chr --distance square 1-ibs flat-missing 


cp /media/elgon/Jesse/HerbariumSeq/plink.mdist ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_allrange.mdist
cp /media/elgon/Jesse/HerbariumSeq/plink.mdist.id ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_allrange.mdist.id




############### MDS for Globe
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange_pruned --allow-extra-chr --cluster --mds-plot 2 --out ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/globalMDS_2



########################### #admixture ##########################

cd /media/elgon/Jesse/HerbariumSeq/

~/Desktop/dist/admixture_linux-1.3.0/admixture stem_outLowMiss_allrange_pruned.bed 2


for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do ~/Desktop/dist/admixture_linux-1.3.0/admixture --cv stem_outLowMiss_allrange_pruned.bed $K | tee log_LowMiss_allrange_pruned${K}.out; done

grep -h CV log_LowMiss_allrange_pruned*.out

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_allrange_pruned.*Q ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/





#### 4. get SNPs for eurasia

cd /media/elgon/Jesse/HerbariumSeq/



##low missing for admixture eurasia
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --geno 0.05 --hwe 1e-50 'midp' --allow-no-sex --mind 0.95 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_eurasianFocus.txt --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia 

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia.imiss ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/




## LD thin eurasia low miss
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia --indep-pairwise 50 10 0.1
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia --extract /media/elgon/Jesse/HerbariumSeq/plink.prune.in --make-bed --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia_pruned




#MDS for Eurasia
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia_pruned --allow-extra-chr --cluster --mds-plot 2 --out ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/EurasiaMDS



#### admixture eurasia ###

~/Desktop/dist/admixture_linux-1.3.0/admixture stem_outLowMiss_Eurasia_pruned.bed 2

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do ~/Desktop/dist/admixture_linux-1.3.0/admixture --cv stem_outLowMiss_Eurasia_pruned.bed $K | tee log_LowMiss_Eurasia_pruned${K}.out; done

grep -h CV log_LowMiss_Eurasia_pruned*.out

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_Eurasia_pruned.*Q ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/




########### 5. get SNPs for S & E Africa

cd /media/elgon/Jesse/HerbariumSeq/


##low missing for admixture eurasia
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --geno 0.05 --hwe 1e-50 'midp' --allow-no-sex --mind 0.95 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_seAfricanFocus.txt --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica.imiss ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/



## LD thin africa low miss
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica --indep-pairwise 50 10 0.1
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica --extract /media/elgon/Jesse/HerbariumSeq/plink.prune.in --make-bed --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica_pruned



####MDS for S & E Africa
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica_pruned --allow-extra-chr --cluster --mds-plot 2 --out ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/seAfricaMDS


#genetic distance thinned
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica_pruned --allow-extra-chr --distance square 1-ibs flat-missing 


cp /media/elgon/Jesse/HerbariumSeq/plink.mdist ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_seAfrica.mdist
cp /media/elgon/Jesse/HerbariumSeq/plink.mdist.id ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_seAfrica.mdist.id



~/Desktop/dist/admixture_linux-1.3.0/admixture stem_outLowMiss_seAfrica_pruned.bed 2


for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; \
do ~/Desktop/dist/admixture_linux-1.3.0/admixture --cv stem_outLowMiss_seAfrica_pruned.bed $K | tee log_LowMiss_seAfrica_pruned${K}.out; done

grep -h CV log_LowMiss_seAfrica_pruned*.out

cp /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_seAfrica_pruned.*Q ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/





##### 6. get GWAS SNPs (not such stringent filtering for missingnesss)

cd /media/elgon/Jesse/HerbariumSeq/


##########

#eurasia
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --allow-no-sex --maf 0.05  --geno 0.2 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_eurasianFocus.txt --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS

#take the stem_out_gwas and make kinship mat
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS --allow-extra-chr --make-rel square 

cp /media/elgon/Jesse/HerbariumSeq/plink.rel ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/
cp /media/elgon/Jesse/HerbariumSeq/plink.rel.id ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/

#GWAS germany
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --allow-no-sex --maf 0.05  --geno 0.2 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forGermanyFocus.txt --pheno /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_pheno.txt --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany

~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany --allow-extra-chr --make-rel square --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany


#GWAS Iberia
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --allow-no-sex --maf 0.05  --geno 0.2 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forIberiaFocus.txt --pheno /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_pheno.txt --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia

~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia --allow-extra-chr --make-rel square  --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia


#GWAS Norway
~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out --allow-extra-chr  --allow-no-sex --maf 0.05  --geno 0.2 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forNorwayFocus.txt --pheno /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_pheno.txt --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway

~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway --allow-extra-chr --make-rel square  --out /media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway




##LD thinning
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss --indep-pairwise 50 10 0.1
#(output indicating number of SNPs targeted for inclusion/exclusion)
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss --extract /media/elgon/Jesse/HerbariumSeq/plink.prune.in --make-bed --out /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_pruned


~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_pruned --missing

cp /media/elgon/Jesse/HerbariumSeq/plink.imiss ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/





#genetic distance thinned
~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_outLowMiss_pruned --allow-extra-chr --distance square 1-ibs flat-missing 


cp /media/elgon/Jesse/HerbariumSeq/plink.mdist* ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/







########################## kinship ##########################

cd /media/elgon/Jesse/HerbariumSeq/

~/Desktop/plink --tfile /media/elgon/Lua/aDNA_2020/FINAL_FILES/stem_out --allow-extra-chr  --geno 0.1 --hwe 1e-50 'midp' --allow-no-sex --mind 0.95 --missing --freq --make-bed --remove  ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps2.txt --out /media/elgon/Jesse/HerbariumSeq/stem_out0.1Miss_forK 


~/Desktop/plink --bfile /media/elgon/Jesse/HerbariumSeq/stem_out0.1Miss_forK --allow-extra-chr --distance square ibs flat-missing 

mv /media/elgon/Jesse/HerbariumSeq/plink.mibs ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink.mibs_stem_out0.1Miss_knship.txt

cp  /media/elgon/Jesse/HerbariumSeq/stem_out0.1Miss_forK.fam ~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/




########################## ##########################




