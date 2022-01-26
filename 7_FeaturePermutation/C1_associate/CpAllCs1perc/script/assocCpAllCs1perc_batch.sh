#!/bin/sh

# seq FIRST STEP LAST
replaceArr=($(seq 825 1 900))

# If files should be named sequentially
#len=${#replaceArr[@]}
#fileArr=($(seq 1 1 ${len}))
fileArr=($(seq 825 1 900))

# Source file
src="assocCpAllCs1perc_template"
#src="job"

# Source file ext.
srcExt=".R"

# Output file
out="assocCpAllCs1perc"

# String to be replaced with var in varArray
orig=FOIREPLACE
lineNum=31 #28
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
