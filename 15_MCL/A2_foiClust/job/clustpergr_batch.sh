#!/bin/sh

#declare -a replaceArr=(0.5 0.6 0.7 0.8 0.9 1)
replaceArr=($(seq 2 1 10))
fileArr=($(seq 2 1 10))

# If files should be named sequentially
#len=${#replaceArr[@]}
#fileArr=($(seq 1 1 ${len}))

# seq FIRST STEP LAST
#varArray=($(seq 1 1 21))
#varArray+=("X")

# Source file
#src="A1_shuffling"
#src="m2c${chr}"
src="clustpergr"
#src="job"

# Source file ext.
srcExt=".sh"

# Output file
#out="m2c"
#out="m2c${chr}cp"
out="clustpergr"

# String to be replaced with var in varArray
orig=REPLACE
lineNum=18
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

# Submitting many jobs

#done
varArray=($(seq 2 1 10))

#declare -a varArray=("Co" "Hi" "Lu" "LV" "RV" "Ao" "PM" "Pa" "Sp" "Li" "SB" "AG"
#"Ov" "Bl" "MesC" "MSC" "NPC" "TLC" "ESC" "LC" "FC")

for var in "${varArray[@]}"
do

qsub clustpergr${var}.sh

done


