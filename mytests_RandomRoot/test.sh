#!/bin/bash

in='../Dataset/37-taxon/noscale.25g.500b/R1/all_gt.tre'
true_tre='../Dataset/37-taxon/true_tree_trimmed'
st_in='all_gt_randroot.tre'
st_out='stelar_output.tre' 


python3 randroot.py $in $st_in
java -jar ../scripts/stelar.jar -i $st_in -o $st_out

index_to_extract=2
tuple=$(python3 ../RF/getFpFn.py -e $st_out -t $true_tre)
RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
echo $RFdistance