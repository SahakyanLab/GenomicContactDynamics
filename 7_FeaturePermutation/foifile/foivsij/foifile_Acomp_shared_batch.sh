#!/bin/sh

#declare -a replaceArr=(0.5 0.6 0.7 0.8 0.9 1)
replaceArr=("Co" "Hi" "Lu" "LV" "RV" "Ao" "PM" "Pa" "Sp" "Li" "SB" "AG"
"Ov" "Bl" "MesC" "MSC" "NPC" "TLC" "ESC" "LC" "FC")
len=${#replaceArr[@]}
fileArr=($(seq 154 1 174)) #${len}))

# Source file
src="foifile_Acomp_shared_template"
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
