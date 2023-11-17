#!/bin/bash
# set -ex
# Local files

# Outer folders with taxa numbers
folders=(  11-taxon 15-taxon 37-taxon 48-taxon )

# Replicates
R=20

# Diffrerent inner folders for each taxa number
innerFolderNames11=(estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS )
innerFolderNames37=(noscale.25g.500b noscale.50g.500b noscale.100g.500b noscale.200g.250b noscale.200g.500b noscale.200g.1000b noscale.200g.1500b noscale.200g.true noscale.400g.500b noscale.800g.500b scale2d.200g.500b scale2u.200g.500b )
innerFolderNames15=(100gene-100bp 100gene-1000bp 100gene-true 1000gene-100bp 1000gene-1000bp 1000gene-true)
innerFolderNames48=(0.5X-1000-500 1X-25-500 1X-50-500 1X-100-500 1X-200-500 1X-500-500 1X-1000-500 2X-1000-500)




declare -A sum_diffs
declare -A sum_rfs
declare -A count_diffs

mkdir -p ../Dataset/Summary_Astral
timediffcsv_main=../Dataset/Summary_Astral/average_diffs.csv
rf_csv_main=../Dataset/Summary_Astral/average_rfs.csv
echo "Folder,AverageDiff" > $timediffcsv_main
printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_diff" > $timediffcsv_main
echo "Folder,Method,AverageRF" > $rf_csv_main
printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv_main

for folder in ${folders[@]}
do
	echo $folder
    mkdir -p ../Dataset/$folder/Summary
	timediffcsv=../Dataset/$folder/Summary/average_diffs_astral.csv
	rf_csv=../Dataset/$folder/Summary/average_rfs_astral.csv
	echo "Folder,AverageDiff" > $timediffcsv
	printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_diff" > $timediffcsv
	echo "Folder,Method,AverageRF" > $rf_csv
	printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv
	
	


	if [ $folder == "11-taxon" ];then
		innerFolderNames=("${innerFolderNames11[@]}")
		R=20
	
	elif [ $folder == "15-taxon" ];then
		innerFolderNames=("${innerFolderNames15[@]}")
		R=10
	
	elif [ $folder == "37-taxon" ];then
		innerFolderNames=("${innerFolderNames37[@]}")
		R=20
	elif [ $folder == "48-taxon" ];then
		innerFolderNames=("${innerFolderNames48[@]}")
		R=20
	fi
	for inner_folder in ${innerFolderNames[@]}
	do
    sum_diffs=0
	count_diffs=0
	sum_rfs=0
    for ((j=1;j<=$R;j++))
		do
			 

			    mkdir -p ../Dataset/$folder/$inner_folder/R$j/astral_outputs/


				input=../Dataset/$folder/$inner_folder/R$j/all_gt.tre
				output=../Dataset/$folder/$inner_folder/R$j/astral_outputs/astral_output.tre
			
				true_tree=../Dataset/$folder/true_tree_trimmed
				touch $output
				echo $input
			
			
				START=$(date +%s.%N)			 
				java -jar astral.5.7.8.jar -i $input -o $output
				END=$(date +%s.%N)
				DIFF=$(echo "$END - $START" | bc)
				echo $DIFF

				index_to_extract=2
				tuple=$(python3 ../RF/getFpFn.py -e $output -t $true_tree)
				RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
				echo $RFdistance

				sum_diffs=$(echo "${sum_diffs} + $DIFF" | bc)
				sum_rfs=$(echo "${sum_rfs} + $RFdistance" | bc)
                ((count_diffs++))
				# break

            
			
		
			 rm -r ../Dataset/$folder/$inner_folder/R$j/astral_outputs
			# break

		done 

        if [ ${count_diffs} -gt 0 ]; then
                average_diff=$(echo "${sum_diffs} / ${count_diffs}" | bc -l)
                printf "Average DIFF for astral: %f\n"  "$average_diff"
                # Uncomment the following line if you want to write the average to a file
                printf "astral,%s,%s,,%.6f\n" "$folder" "$inner_folder"  "$average_diff" >> $timediffcsv
				printf "astral,%s,%s,,%.6f\n" "$folder" "$inner_folder"  "$average_diff" >> $timediffcsv_main
                average_rf=$(echo "${sum_rfs} / ${count_diffs}" | bc -l)
				printf "Average RF for astral: %f\n"  "$average_rf"
                
                printf "astral,%s,%s,,%.6f\n" "$folder" "$inner_folder"  "$average_rf" >> $rf_csv
				printf "astral,%s,%s,,%.6f\n" "$folder" "$inner_folder"  "$average_rf" >> $rf_csv_main
              

        fi
        printf "\n" >> $timediffcsv
		printf "\n" >> $timediffcsv_main
		
		printf "\n" >> $rf_csv
		printf "\n" >> $rf_csv_main

        done

        done