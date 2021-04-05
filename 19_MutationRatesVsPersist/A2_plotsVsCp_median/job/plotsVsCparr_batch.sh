#!/bin/sh

replaceArr=("Tmut" "Nmsite" "TmutDIVNmsite" "Nmsitenorm" "numWTSEQ")

# If files should be named sequentially
len=${#replaceArr[@]}
fileArr=($(seq 1 1 ${len}))

# Source file
src="plotsVsCparr"

# Source file ext.
srcExt=".sh"

# Output file
out="calc"

# String to be replaced with var in varArray
orig=CALCREPLACE
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
done
###################################################################################


