#!/bin/sh

#array1=("numWTSEQ" "Tmut" "Tmutnorm" "Nmsite" "Nmsitenorm" "TmutDIVNmsite")
array1=("Tmutnorm")
prefix1="calc"

array2=($(seq 1 1 1400)) #arr2.repl
prefix2="comb"
len2=${#array2[@]}
fileArr2=($(seq 1 1 ${len2}))

srcfile="/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist/A2_plotdataVsCp_median/script/A2_plotdataVsCp_template.R"
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
for ar1 in "${array1[@]}"
do

outNme=${prefix1}${ar1}
cp ${srcfile} "${outNme}.R"

ar1=${ar1}
sed -e "s/arr1.repl/${ar1}/" "${outNme}.R" > temp3759157105${ar1}
mv temp3759157105${ar1} "${outNme}.R"

#for ar2 in "${array2[@]}"
for (( i=0; i<${len2}; i++ ))
do

fileind2=${fileArr2[$i]}
ar2=${array2[$i]}

outNme1=${outNme}${prefix2}${fileind2}
# Make a copy of the source file and rename to output filename
cp "${outNme}.R" "${outNme1}.R"
# Replace line
sed -e "s/arr2.repl/${ar2}/" "${outNme1}.R" > temp3759157105${ar1}${ar2}
mv temp3759157105${ar1}${ar2} "${outNme1}.R"
done # array2

done
   
