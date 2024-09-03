#!/bin/bash
for f in *.bedgraph
do
 echo "Processing $f"
 grep -E 'chr19|chr7|chr11' $f > ${f%.bedgraph}.chr-spec.bedgraph
done
