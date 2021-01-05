#!/bin/sh

#declare -a replaceArr=(0.5 0.6 0.7 0.8 0.9 1)
replaceArr=("ESC" "FC" "LC")
len=${#replaceArr[@]}
fileArr=($(seq 100 1 102))

# Source file
src="foifile_H3K9me3_template"
# Source file ext.
srcExt=""

# Output file
out="foifile"

# String to be replaced with var in varArray
orig=REPLACE
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
