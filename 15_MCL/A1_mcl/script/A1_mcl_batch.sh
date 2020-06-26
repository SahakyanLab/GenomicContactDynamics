#!/bin/sh

# seq FIRST STEP LAST
# Chr array
#array1=($(seq 1 1 23))
declare -a array1=("Cs" "Cp")
# Cp array
#array2=($(seq 1 1 21))
declare -a array2=(6 7 8 9)
srcfile="/Users/ltamon/DPhil/GenomicContactDynamics/21_MCL/A1_mcl/script/A1_mcl.R"
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
for ar1 in "${array1[@]}"
do

outNme="m2w${ar1}"
cp ${srcfile} "${outNme}.R"

# if ar1=23, replace with 'X'
#if [ "$ar1" -eq "23" ]
#then
#ar1="X";
#else
#ar1=$ar1;
#fi

sed -e "s/arr1.repl/${ar1}/" "${outNme}.R" > temp3759157105${ar1}
mv temp3759157105${ar1} "${outNme}.R"

for ar2 in "${array2[@]}"
do

outNme1="${outNme}gr${ar2}"
# Make a copy of the source file and rename to output filename
cp "${outNme}.R" "${outNme1}.R"
# Replace line
sed -e "s/arr2.repl/${ar2}/" "${outNme1}.R" > temp3759157105${ar1}${ar2}
mv temp3759157105${ar1}${ar2} "${outNme1}.R"
done # array2

done

