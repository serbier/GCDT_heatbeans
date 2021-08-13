#!/bin/bash
WD=./../..
NGSEP=/home/dariza/software/NGSEP/NGSEPcore_4.0.0.jar

#filter vcf to nonmiss 
file=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
VCF=${WD}/data/genotype/${file}.vcf.gz
OUT=${WD}/SNP_GWAS/genotype/${file}

#bcftools view -i 'MAF<=0.03' ${VCF}| bgzip > ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3.vcf.gz
#echo "SNP Number:" `bcftools view -H ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3.vcf.gz| wc -l`

#Hapmap transformation
java -jar ${NGSEP} VCFConverter -i ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3.vcf.gz -o ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3 -hapmap


BLUPs=${WD}/Phenotipyc_analysis/GCDT_MERGED/results/merged_GCDT_blups_Tepary_GWAS.csv
BLUEs=${WD}/Phenotipyc_analysis/GCDT_MERGED/results/merged_GCDT_blues_Tepary_GWAS.csv

mkdir -p ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUPs

nohup Rscript ${WD}/SNP_GWAS/src/GAPIT_GWAS_MLM_BLINK.R -geno ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3_hmp.txt\
  -pheno ${BLUPs}\
  -pca 5\
  -modelSelection TRUE\
  -out ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUPs > ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUPs/log_file_$(date '+%Y-%m-%d').out &

mkdir -p ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUEs

nohup Rscript ${WD}/SNP_GWAS/src/GAPIT_GWAS_MLM_BLINK.R -geno ${WD}/SNP_GWAS/genotype/GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_MAF_less_3_hmp.txt\
  -pheno ${BLUEs}\
  -pca 5\
  -modelSelection TRUE\
  -out ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUEs > ${WD}/SNP_GWAS/results/MAF_less_3_GAPIT_PCA/BLUEs/log_file_$(date '+%Y-%m-%d').out &
