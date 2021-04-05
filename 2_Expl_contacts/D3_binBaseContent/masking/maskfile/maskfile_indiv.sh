#!/bin/sh

lines=`cat maskfile_master`
outid="maskfile"

i=0

for line in $lines; do
((i++))
echo "$line" > ${outid}${i}
done
