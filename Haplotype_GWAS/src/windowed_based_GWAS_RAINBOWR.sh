#!/bin/bash
#In this script we perform a SNP-set based GWAS with R package RAINBOWR
#RAINBOWR(Reliable Association INference By Optimizing Weights with R)
WD=./../..
fileName=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
vcf=${WD}/data/genotype/${fileName}.vcf.gz

plink --vcf ${vcf} --allow-extra-chr --double-id --make-bed  --out ${WD}/Haplotype_GWAS/genotype/${fileName}

