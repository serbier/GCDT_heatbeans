
#!/bin/bash
WD=./
file=GCDT_plates_fix_Vulgaris_v2_Bi_Q40_Dp3_imiss83_NS99_MAF2_He2_annotated
NGSEP=${WD}/bin/NGSEPcore_4.0.0.jar

java  -Xmx10G -jar ${NGSEP} VCFIntrogressionAnalysis \
	-i ${WD}/genotype/${file}.vcf.gz\
	-p ${WD}/data/NGSEP/groups.txt\
	-o GCDT_Acu_vul\
	-g 80\
	-m 0.4\
	-w 50\
	-v 0\
	-s 30\
	-c\
	-u

mv GCDT_Acu_vul* ${WD}/results/NGSEP/