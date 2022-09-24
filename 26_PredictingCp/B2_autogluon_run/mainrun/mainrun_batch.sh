#!/bin/sh

array1=("['chr' + str(chr_num) for chr_num in [1]]") #("['chr' + str(chr_num) for chr_num in list(range(18))]") #arr1.repl
len1=${#array1[@]}
fileArr1=($(seq 1 1 ${len1}))
prefix1="chrid"

array2=("grp.compl.|grp.kmer3." "grp.compl." "grp.kmer3." "grp.compl.|grp.kmer3.|grp.kmer1.|grp.GC.")  #arr2.repl
len2=${#array2[@]}
fileArr2=($(seq 1 1 ${len2}))
prefix2="feat"

#srcfile='/Users/ltamon/SahakyanLab/GenomicContactDynamics/26_PredictingCp/B2_autogluon_run/mainrun/mainrun_template.py'
srcfile='/project/sahakyanlab/ltamon/SahakyanLab/GenomicContactDynamics/26_PredictingCp/B2_autogluon_run/mainrun/mainrun_template.py'
srcExt='.py'
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
for (( i=0; i<${len1}; i++ ))
do
	
	fle1=${fileArr1[$i]}
	outNme1=${prefix1}${fle1}
	cp ${srcfile} "${outNme1}${srcExt}"
	
	ar1=${array1[$i]}
	sed -e "s/arr1.repl/${ar1}/" "${outNme1}${srcExt}" > temp3759157105${fle1}
	mv temp3759157105${fle1} "${outNme1}${srcExt}"
	
	for (( j=0; j<${len2}; j++ ))
	do
		
		fle2=${fileArr2[$j]}
		outNme2=${outNme1}${prefix2}${fle2}
		cp "${outNme1}${srcExt}" "${outNme2}${srcExt}"
		
		ar2=${array2[$j]}
		sed -e "s/arr2.repl/${ar2}/" "${outNme2}${srcExt}" > temp3759157105${fle1}${fle2}
		mv temp3759157105${fle1}${fle2} "${outNme2}${srcExt}"
		
	done # array2
	
done # array1

