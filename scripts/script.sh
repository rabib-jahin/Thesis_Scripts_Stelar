#!/bin/bash
# set -ex

timeFile="11-taxa-time.csv" # will keep the records of running time
RFRateFile="11-RFRateFile.csv" # will keep the records of RF rates


folders=( 11-taxon  )
innerFolderNames=(estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS )
methodNames=( MP MV )
# Headers of the csv files
#echo "ilslevel, genetrees, basepair, methodName, modelConditionName,Robinson-Foulds distance" > $RFRateFile
#echo "name,time" > $timeFile

for folder in ${folders[@]}
do
	echo $folder
	speciesTree=../Dataset/$folder/true_tree_trimmed
	echo $speciesTree
	
	for inner_folder in ${innerFolderNames[@]}
	do
	  
	 
		for j in {1..20}
		do
			 
			gt_folder=$inner_folder/R$j
			
			gt=$gt_folder/all_gt.tre
			

			# if [[ -f "$speciesTree.stripped" ]];
			# then 
			# 	continue
			# fi
			for method in  ${methodNames[@]}

			do
 
			                mkdir -p ../Dataset/$folder/$gt_folder/stelar_inputs
					file=../Dataset/$folder/$gt
					START=$(date +%s.%N)
					python3 utils.py $file temp.tre
					python3 ../fastroot/MinVar-Rooting/FastRoot.py -m $method -i temp.tre  -o ../Dataset/$folder/$gt_folder/stelar_inputs/stelar_input_$method.tre 
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
					# cp gt.best.of.10.tre  $speciesTree
					# rm -r gt*
					echo $DIFF
					rm temp.tre
				

             done
			# ./strip_edge_support2.pl -i $speciesTree -o $speciesTree.stripped
			# RFdistance=$(python2 getFpFn.py -e $speciesTree.stripped -t true.tree.strip | sed 's/.//; s/,//' |awk '{print $3}')

			# printf $folder | tr '.' ',' >> $RFRateFile; echo ,$method,R$j,"${RFdistance//)}" >> $RFRateFile
			# echo $folder,R$j,$method,$DIFF >> $timeFile
			break
		done
		break
	done
	break
done

