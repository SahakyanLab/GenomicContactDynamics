#!/bin/sh

replaceArr=($(seq 9 1 12))
fileArr=($(seq 9 1 12))

# Source file
src="C2_generate_map_template"

# Source file ext.
srcExt=".R"

# Output file
out="genmap"

# String to be replaced with var in varArray
orig=PARAMREPLACE
lineNum=37
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
