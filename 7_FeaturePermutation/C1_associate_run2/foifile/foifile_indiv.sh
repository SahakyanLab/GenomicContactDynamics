#!/bin/sh

lines=`cat foifile_agerank`
outid="foifile"

i=59

for line in $lines; do
((i++))
echo "$line" > ${outid}${i}
done
