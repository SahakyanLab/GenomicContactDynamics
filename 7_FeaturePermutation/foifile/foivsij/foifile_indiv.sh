#!/bin/sh

lines=`cat foifile24`
outid="foifile"
fileArr=($(seq 27 1 32))

i=0

for line in $lines; do
fileind=${fileArr[$i]}
echo "$line" > ${outid}${fileind}
((i++))
done
