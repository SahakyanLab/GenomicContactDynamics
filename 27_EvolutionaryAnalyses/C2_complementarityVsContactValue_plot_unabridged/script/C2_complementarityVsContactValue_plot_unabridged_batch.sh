#!/bin/sh

replaceArr=(median.consCp max.consCp min.consCp mean.consCp MODE.consCp)
fileArr=($(seq 1 1 5))
src="C2_complementarityVsContactValue_plot_unabridged_template"
srcExt=".R"
out="UAplot"
orig=CPREPLACE
lineNum=33
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
len=${#fileArr[@]}

for (( i=0; i<${len}; i++ ))
do

fileind=${fileArr[$i]}
var=${replaceArr[$i]}

# Full filename of output file
outNme="${out}${fileind}${srcExt}"
# Make a copy of the source file and rename to output filename
cp ${src}${srcExt} ${outNme}
# Replace line
sed -e "${lineNum}s/${orig}/${var}/" ${outNme} > tempfile
mv tempfile ${outNme}
done # varArray
###################################################################################

