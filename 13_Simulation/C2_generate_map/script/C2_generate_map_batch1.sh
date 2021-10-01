#!/bin/sh

# declare -a replaceArr=(0.5 0.6 0.7 0.8 0.9 1)
# replaceArr=(Co Hi Lu LV RV Ao PM Pa Sp Li SB AG Ov Bl MesC MSC NPC TLC ESC FC LC)
replaceArr=("Cs.raw" "Cs.norm" "Cp" "CII.disc.kmer.5" "CII.cont.kmer.5")
# If files should be named sequentially
len=${#replaceArr[@]}
fileArr=($(seq 1 1 ${len}))

# Source file
src="C2_generate_map_template"

# Source file ext.
srcExt=".R"

# Output file
out="LCmetric"

# String to be replaced with var in varArray
orig=METRICREPLACE
lineNum=52
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

#declare -a varArray=("Co" "Hi" "Lu" "LV" "RV" "Ao" "PM" "Pa" "Sp" "Li" "SB" "AG"
#"Ov" "Bl" "MesC" "MSC" "NPC" "TLC" "ESC" "LC" "FC")

#for var in "${varArray[@]}"
#do

#sbatch ${var}genmaparr.sh

#done



