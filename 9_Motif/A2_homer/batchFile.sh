#!/bin/sh

declare -a replaceArr=(0.5 0.6 0.7 0.8 0.9 1)

# seq FIRST STEP LAST
#varArray=($(seq 1 1 21))
#varArray+=("X")

# Source file
#src="A1_shuffling"
#src="m2c${chr}"
src="A2_homer"
#src="job"

# Source file ext.
srcExt=".R"

# Output file
#out="m2c"
#out="m2c${chr}cp"
out="A2_homer_bgrCp1_"

# String to be replaced with var in varArray
orig=REPLACE
lineNum=40
###################################################################################
# MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE # # MAIN CODE #
###################################################################################
len=${#replaceArr[@]}

fileArr=($(seq 1 1 ${len}))

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

