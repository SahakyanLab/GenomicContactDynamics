#!/bin/sh

lines=`cat foifile_list`
outid="foifile"
fileArr=($(seq 1 1 11))

i=0

for line in $lines; do
fileind=${fileArr[$i]}
echo "$line" > ${outid}${fileind}
((i++))
done
