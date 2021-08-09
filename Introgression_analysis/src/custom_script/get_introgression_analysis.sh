#!/bin/bash
WD=./../..
file=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
pop_parents_vul=${WD}/data/custom_script/pop_parents_vul.txt
pop_parents_acut=${WD}/data/custom_script/pop_parents_acut.txt

#get frequencies for population pop_parents_vul
#vcftools --gzvcf ${WD}/genotype/${file}.vcf.gz --freq --keep ${pop_parents_vul} --out ${WD}/data/custom_script/${file}_pop_parents_vul
#get frequencies for population pop_parents_acut
#vcftools --gzvcf ${WD}/genotype/${file}.vcf.gz --freq --keep ${pop_parents_acut} --out ${WD}/data/custom_script/${file}_pop_parents_acut

#get chrom pos of vcf file for genetic position 
#bcftools query -f '%CHROM\t%POS\n' ${WD}/genotype/${file}.vcf.gz > ${WD}/data/custom_script/${file}_positions.txt
#Rscript ${WD}/src/custom_script/getGenPos.R 

#get TGT vcf 
#bcftools query -f '%CHROM\t%POS\t[%TGT\t]\n' ${WD}/genotype/${file}.vcf.gz > ${WD}/data/custom_script/${file}_TGT.vcf
#get samples of vcf
#bcftools query -l ${WD}/genotype/${file}.vcf.gz > ${WD}/data/custom_script/samples.txt


function printThreads {

  tmp=$1
  lst=(`cat ${tmp} | tr '\n' ' '`)

  myNum=0
    echo -e 'File '${tmp}' contains the following samples:\n'${lst[@]}
    myNum=`expr ${myNum} + ${#lst[@]}`
    echo -e 'No. of samples assigned: '${myNum}'\n'
}


function assignThreads {

  # 1st argument is a bash array with sample names
  arr=("$@")
  # 2nd argument is the number of threads to be used in the machine
  nt=${arr[-1]}
  unset 'arr[${#arr[@]}-1]'

  for i in ${!arr[@]}
  do

    spl=`expr ${i} % ${nt}`
    echo ${arr[${i}]} >> ${WD}/src/custom_script/tmpDir/tmpList_${spl}.tmp
  
  done
}

list=(`cat ${WD}/data/custom_script/samples.txt | tr '\n' ' '`)
mkdir tmpDir
rm ${WD}/src/custom_script/tmpDir/*.tmp
nThreads=5

echo -e '\nStarting introgression mapping in '${WD}' '$(date)''

rm ${WD}/src/custom_script/tmpDir/*.tmp


assignThreads ${list[@]} ${nThreads}

for tmpF in ${WD}/src/custom_script/tmpDir/tmpList*.tmp
do
    tList=(`cat ${tmpF} | tr '\n' ' '`)
    printThreads ${tmpF}

  ( for sample in ${tList[@]}
    do sleep 1
      echo $sample
      python introgression_profile.py --s $sample --vcf ${WD}/data/custom_script/${file}_TGT.vcf --genpos ${WD}/data/custom_script/${file}_gen_pos.csv --out ${WD}/results/custom_script/outs/
    done ) &

done
wait