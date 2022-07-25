#!/bin/sh

replaceArr=($(seq 1 1 21))
fileArr=($(seq 1 1 21))

# Source file
src="ijconsensus_template"

# Source file ext.
srcExt=".sh"

# Output file
out="ijcons"

# String to be replaced with var in varArray
orig=CPREPLACE
lineNum=16
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

