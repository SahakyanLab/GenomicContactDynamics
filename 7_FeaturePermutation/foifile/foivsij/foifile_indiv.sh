#!/bin/sh

lines=`cat foifile_isochore`
outid="foifile"
fileArr=($(seq 177 1 188))

i=0

for line in $lines; do
fileind=${fileArr[$i]}
echo "$line" > ${outid}${fileind}
((i++))
done
