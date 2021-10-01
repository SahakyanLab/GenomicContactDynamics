#!/bin/sh

#array1=(Co Hi Lu LV RV Ao PM Pa Sp Li SB AG Ov Bl MesC MSC NPC TLC ESC FC LC)
array1=(AG Ao Bl Co ESC FC Hi LC Li Lu LV MesC MSC NPC Ov Pa PM RV SB Sp TLC)
array2=(SIM.non.set3.2.20.5.1.0 SIM.non.set3.2.20.5.1.5 SIM.int.set3.2.20.5.1.0 SIM.int.set3.2.20.5.1.5 SIM.non.set3.2.10.5.1.0 SIM.non.set3.2.10.5.1.5 SIM.int.set3.2.10.5.1.0 SIM.int.set3.2.10.5.1.5)
array3=(Cs.norm)
arrID=(ct subj ref)
srcfile="/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation/B1_compare/Cs.norm/compare/script/compare_template.R"
repStr=(CTREPLACE SUBJREPLACE REFREPLACE)
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
arr1len=${#array1[@]}
arr2len=${#array2[@]}
arr3len=${#array3[@]}

for (( i=0; i<${arr1len}; i++ ))
do

a1=${array1[${i}]}
fInd1=${i}
outNme1="${arrID[0]}$((++fInd1))"
cp ${srcfile} "${outNme1}.R"
sed -e "s/${repStr[0]}/${a1}/" "${outNme1}.R" > temp341874${i}
mv temp341874${i} "${outNme1}.R"

for (( j=0; j<${arr2len}; j++ ))
do

a2=${array2[${j}]}
fInd2=${j}
outNme2="${outNme1}${arrID[1]}$((++fInd2))"
cp "${outNme1}.R" "${outNme2}.R"
sed -e "s/${repStr[1]}/${a2}/" "${outNme2}.R" > temp341874${i}${j}
mv temp341874${i}${j} "${outNme2}.R"

for (( k=0; k<${arr3len}; k++ ))
do

a3=${array3[${k}]}
fInd3=${k}
outNme3="${outNme2}${arrID[2]}$((++fInd3))"
cp "${outNme2}.R" "${outNme3}.R"
sed -e "s/${repStr[2]}/${a3}/" "${outNme3}.R" > temp341874${i}${j}${k}
mv temp341874${i}${j}${k} "${outNme3}.R"

done # array3

rm "${outNme2}.R"

done # array2

rm "${outNme1}.R"

done # array1

