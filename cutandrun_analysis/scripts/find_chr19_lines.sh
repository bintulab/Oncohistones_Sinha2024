#!/bin/bash
for f in *_dedup.scaled.bedgraph
do
 echo "Processing $f"
 grep 'chr19' $f > ${f%_edited_dedup.scaled.bedgraph}_chr19_dedup.scaled.bedgraph
done
