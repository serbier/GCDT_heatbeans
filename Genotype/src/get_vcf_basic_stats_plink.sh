#!/bin/bash
WD=./../..
fileName=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
vcf=${WD}/data/genotype/${fileName}.vcf.gz

#convert to plink standard format
plink --vcf ${vcf} --allow-extra-chr --double-id --recode --out ${WD}/data/processed_data/${fileName}

#get missingess stats
plink --file ${WD}/data/processed_data/${fileName} --missing --nonfounders --allow-extra-chr --out ${WD}/data/results/missingness

#get allele frequencies
plink --file ${WD}/data/processed_data/${fileName} --freq --nonfounders --allow-extra-chr --out ${WD}/data/results/allele_freqs


