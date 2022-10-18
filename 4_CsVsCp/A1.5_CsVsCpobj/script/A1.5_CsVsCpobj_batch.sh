#!/bin/sh

replaceArr=(Co Hi Lu LV RV Ao PM Pa Sp Li SB AG Ov Bl MesC MSC NPC TLC ESC FC LC)
lenrep=${#replaceArr[@]}
fileArr=($(seq 1 1 ${lenrep}))

# Source file
src="A1.5_CsVsCpobj_template"
srcExt=".R"

# Output file
out="objct"

# String to be replaced with var in varArray
orig=REPLACE
lineNum=31
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
