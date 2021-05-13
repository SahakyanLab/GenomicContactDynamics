#!/bin/sh

replaceArr=("All" "C>A" "C>G" "C>T" "T>A" "T>C" "T>G")

# If files should be named sequentially
len=${#replaceArr[@]}
fileArr=($(seq 1 1 ${len}))

# Source file
src="A1_mutCalcPerBin_template"

# Source file ext.
srcExt=".R"

# Output file
out="mutCalcPerBin"

# String to be replaced with var in varArray
orig=MUTREPLACE
lineNum=44
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
done
###################################################################################


