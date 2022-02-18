#!/bin/sh

array1=($(seq 1 1 11)) #arr1.repl, note that if array1 has 23 --> X, disable this if
# array1 is not for chr
prefix1="foi"
array2=(kmer align)  #arr2.repl
prefix2="type"
srcfile="/Users/ltamon/DPhil/GCD_polished/12_Shuffling/B3_compare/script/B3_compare_template.R"
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
for ar1 in "${array1[@]}"
do

outNme=${prefix1}${ar1}
cp ${srcfile} "${outNme}.R"

# if ar1=23, replace with 'X'
if [ "$ar1" -eq "23" ]
then
ar1="X";
else
ar1=$ar1;
fi

ar1=${ar1}
sed -e "s/arr1.repl/${ar1}/" "${outNme}.R" > temp3759157105${ar1}
mv temp3759157105${ar1} "${outNme}.R"

for ar2 in "${array2[@]}"
do

outNme1=${outNme}${prefix2}${ar2}
# Make a copy of the source file and rename to output filename
cp "${outNme}.R" "${outNme1}.R"
# Replace line
sed -e "s/arr2.repl/${ar2}/" "${outNme1}.R" > temp3759157105${ar1}${ar2}
mv temp3759157105${ar1}${ar2} "${outNme1}.R"
done # array2

done

