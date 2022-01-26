#!/bin/sh

lines=`cat foifile_raw_TBA`
outid="foifile"

i=824

for line in $lines; do
((i++))
echo "$line" > ${outid}${i}
done
