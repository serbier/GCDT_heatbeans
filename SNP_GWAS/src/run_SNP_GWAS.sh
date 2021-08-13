#!/bin/bash
WD=./../..
NGSEP=/home/dariza/software/NGSEP/NGSEPcore_4.0.0.jar

#Convert vcf to hapmap with NGSEP
file=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
VCF=${WD}/data/genotype/${file}.vcf.gz
OUT=${WD}/SNP_GWAS/genotype/${file}
#java -jar ${NGSEP} VCFConverter -i ${VCF} -o ${OUT} -hapmap

#RUN GWAS Analysis using GAPIT package implementing MLM & BLINK models with PCA method
#for adjust population structure with BLUPs and BLUEs as phenotypic information 

BLUPs=${WD}/Phenotipyc_analysis/GCDT_MERGED/results/merged_GCDT_blups_Tepary_GWAS.csv
BLUEs=${WD}/Phenotipyc_analysis/GCDT_MERGED/results/merged_GCDT_blues_Tepary_GWAS.csv

mkdir -p ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUPs

nohup Rscript ${WD}/SNP_GWAS/src/GAPIT_GWAS_MLM_BLINK.R -geno ${OUT}_hmp.txt\
  -pheno ${BLUPs}\
  -pca 5\
  -modelSelection TRUE\
  -out ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUPs > ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUPs/log_file_$(date '+%Y-%m-%d').out &

mkdir -p ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUEs

nohup Rscript ${WD}/SNP_GWAS/src/GAPIT_GWAS_MLM_BLINK.R -geno ${OUT}_hmp.txt\
  -pheno ${BLUEs}\
  -pca 5\
  -modelSelection TRUE\
  -out ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUEs > ${WD}/SNP_GWAS/results/GAPIT_PCA/BLUEs/log_file_$(date '+%Y-%m-%d').out &
