#!/bin/sh

replaceArr=($(seq 1 1 105))
fileArr=($(seq 1 1 105))

# Source file
src="A3_mut_contact_Cp_fast_template"
srcExt=".R"

# Output file
out="fast"

# String to be replaced with var in varArray
orig=INDREPLACE
lineNum=58
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