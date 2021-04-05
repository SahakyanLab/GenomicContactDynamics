#!/bin/sh

declare -a replaceArr=(1 2 3 4)
len=${#replaceArr[@]}
fileArr=($(seq 1 1 ${len}))

# Source file
src="binBaseContMask_template"
# Source file ext.
srcExt=".sh"

# Output file
out="binBaseContMask"

# String to be replaced with var in varArray
orig=MASKREPLACE
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
sed -e "s/${orig}/${var}/" ${outNme} > tempfile
mv tempfile ${outNme}

done # varArray
###################################################################################
