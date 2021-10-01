#!/bin/sh

lines=`cat foifile_master`
outid="foifile"

i=0

for line in $lines; do
((i++))
echo "$line" > ${outid}${i}
done
