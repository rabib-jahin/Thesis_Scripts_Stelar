#!/bin/bash
# set -ex

timeFile="11-taxa-time.csv" # will keep the records of running time
RFRateFile="11-RFRateFile.csv" # will keep the records of RF rates


folders=( 11-taxon 15-taxon 37-taxon)
R=20
innerFolderNames11=(estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS )
innerFolderNames37=(noscale.25g.500b noscale.50g.500b noscale.100g.500b noscale.200g.250b noscale.200g.500b noscale.200g.1000b noscale.200g.1500b noscale.200g.true noscale.400g.500b noscale.800g.500b scale2d.200g.500b scale2u.200g.500b )
innerFolderNames15=(100gene-100bp 100gene-1000bp 100gene-true 1000gene-100bp 1000gene-1000bp 1000gene-true)
methodNames=(MP)
# Headers of the csv files
#echo "ilslevel, genetrees, basepair, methodName, modelConditionName,Robinson-Foulds distance" > $RFRateFile
#echo "name,time" > $timeFile

for folder in ${folders[@]}
do
    
	echo $folder
	speciesTree=../Dataset/$folder/true_tree_trimmed
	echo $speciesTree
	if [ $folder == "11-taxon" ];then
	innerFolderNames=("${innerFolderNames11[@]}")
	
	elif [ $folder == "15-taxon" ];then
	innerFolderNames=("${innerFolderNames15[@]}")
	R=10
	
	elif [ $folder == "37-taxon" ];then
	innerFolderNames=("${innerFolderNames37[@]}")
	

	
	else
 echo ""
	fi
	
	echo  ${innerFolderNames[@]}
	for inner_folder in ${innerFolderNames[@]}
	do
	
	  
	 
		for ((j=1;j<=$R;j++))
		do
			 
			gt_folder=$inner_folder/R$j
			echo "sss"
			
			gt=$gt_folder/all_gt.tre
			

			# if [[ -f "$speciesTree.stripped" ]];
			# then 
			# 	continue
			# fi
			mkdir -p ../Dataset/$folder/$gt_folder/stelar_inputs
			for method in  ${methodNames[@]}

			do
 
			                
					file=../Dataset/$folder/$gt
					START=$(date +%s.%N)
					echo here
					python3 utils.py $file temp.tre
					python3 ../fastroot/MinVar-Rooting/FastRoot.py -m $method -i temp.tre  -o ../Dataset/$folder/$gt_folder/stelar_inputs/stelar_input_$method.tre 
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
					# cp gt.best.of.10.tre  $speciesTree
					# rm -r gt*
					echo $DIFF
					# rm temp.tre
				

             done
			# ./strip_edge_support2.pl -i $speciesTree -o $speciesTree.stripped
			# RFdistance=$(python2 getFpFn.py -e $speciesTree.stripped -t true.tree.strip | sed 's/.//; s/,//' |awk '{print $3}')

			# printf $folder | tr '.' ',' >> $RFRateFile; echo ,$method,R$j,"${RFdistance//)}" >> $RFRateFile
			# echo $folder,R$j,$method,$DIFF >> $timeFile
			# break
		done
		# break
	done
	# break
done




declare -A sum_diffs
declare -A sum_rfs
declare -A count_diffs


mkdir -p ../Dataset/$folder/Summary
timediffcsv=../Dataset/$folder/Summary/average_diffs.csv
rf_csv=../Dataset/$folder/Summary/average_rfs.csv
echo "Folder,Method,AverageDiff" > $timediffcsv
printf "%s,%s,%s,%s\n\n" "folder" "inner_folder" "method" "average_diff" > $timediffcsv
echo "Folder,Method,AverageRF" > $rf_csv
printf "%s,%s,%s,%s\n\n" "folder" "inner_folder" "method" "average_rf" > $rf_csv


for folder in ${folders[@]}
do
	echo $folder
	speciesTree=../Dataset/$folder/true_tree_trimmed
	echo $speciesTree
	
	for inner_folder in ${innerFolderNames[@]}
	do
		for method in  ${methodNames[@]}
		do
			# Initialize sum and count for each method
			sum_diffs[$method]=0
			count_diffs[$method]=0
			sum_rfs[$method]=0
		done
	  
	 
		for ((j=1;j<=$R;j++))
		do
			 
			# gt_folder=$inner_folder/R$j
			
			# gt=$gt_folder/all_gt.tre
			mkdir -p ../Dataset/$folder/$inner_folder/R$j/stelar_outputs/
			input_folder=$inner_folder/R$j/stelar_inputs
			output_folder=$inner_folder/R$j/stelar_outputs

			# if [[ -f "$speciesTree.stripped" ]];
			# then 
			# 	continue
			# fi
			mkdir -p ../Dataset/$folder/$gt_folder/stelar_outputs
			for method in  ${methodNames[@]}
			do
				# if [ $method == "stelar" ]
				# then
				# 	echo in $method
				# 	java -jar stelar.jar -i  -o < out-species-tree-name >

				# fi
				input=../Dataset/$folder/$inner_folder/R$j/stelar_inputs/stelar_input_$method.tre
				input_sanitized=../Dataset/$folder/$inner_folder/R$j/stelar_inputs/stelar_input_sanitized$method.tre
				output=../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.tre
				stripped_output=../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.stripped.tre
				true_tree=../Dataset/$folder/true_tree_trimmed

				touch $output
				echo $input
				

				python3 helper.py $input $input_sanitized
			
				START=$(date +%s.%N)
			
						
				
				 
				java -jar stelar.jar -i $input_sanitized -o $output


				
				
				END=$(date +%s.%N)
				DIFF=$(echo "$END - $START" | bc)
				echo $DIFF
				index_to_extract=2
				tuple=$(python3 ../RF/getFpFn.py -e $output -t $true_tree)
				RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
				echo $RFdistance


				sum_diffs[$method]=$(echo "${sum_diffs[$method]} + $DIFF" | bc)
				sum_rfs[$method]=$(echo "${sum_rfs[$method]} + $RFdistance" | bc)
                ((count_diffs[$method]++))
				


				#uncomment for RFScore
				# ./strip_edge_support2.pl -i $output -o $stripped_output
				# RFdistance=$(python2 getFpFn.py -e $stripped_output -t $true_tree | sed 's/.//; s/,//' |awk '{print $3}')
				# echo $RFdistance
				

             done
			# ./strip_edge_support2.pl -i $speciesTree -o $speciesTree.stripped
			# RFdistance=$(python2 getFpFn.py -e $speciesTree.stripped -t true.tree.strip | sed 's/.//; s/,//' |awk '{print $3}')

			# printf $folder | tr '.' ',' >> $RFRateFile; echo ,$method,R$j,"${RFdistance//)}" >> $RFRateFile
			# echo $folder,R$j,$method,$DIFF >> $timeFile
			# break
		done

		

		# Print average DIFF for each method

		for method in  ${methodNames[@]}
        do
            if [ ${count_diffs[$method]} -gt 0 ]; then
                average_diff=$(echo "${sum_diffs[$method]} / ${count_diffs[$method]}" | bc -l)
                printf "Average DIFF for %s: %f\n" "$method" "$average_diff"
                # Uncomment the following line if you want to write the average to a file
                printf "%s,%s,%s,%.6f\n\n" "$folder" "$inner_folder" "$method" "$average_diff" >> $timediffcsv
                average_rf=$(echo "${sum_rfs[$method]} / ${count_diffs[$method]}" | bc -l)
				printf "Average RF for %s: %f\n" "$method" "$average_rf"
                
                printf "%s,%s,%s,%.6f\n\n" "$folder" "$inner_folder" "$method" "$average_rf" >> $rf_csv
              

            fi
        done

		# break
	done
	# break
done

