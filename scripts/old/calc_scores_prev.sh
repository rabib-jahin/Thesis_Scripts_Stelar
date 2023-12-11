#!/bin/bash
# set -ex
# Local files

# Outer folders with taxa numbers: 11-taxon, 15-taxon, 37-taxon, 48-taxon
folders=( 15-taxon )
fresh=1
# Replicates
R=20

# Diffrerent inner folders for each taxa number
innerFolderNames11=(estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS )
innerFolderNames37=(noscale.25g.500b noscale.50g.500b noscale.100g.500b noscale.200g.250b noscale.200g.500b noscale.200g.1000b noscale.200g.1500b noscale.200g.true noscale.400g.500b noscale.800g.500b scale2d.200g.500b scale2u.200g.500b )
innerFolderNames15=(100gene-100bp/estimated-genetrees 100gene-1000bp/estimated-genetrees 100gene-true 1000gene-100bp/estimated-genetrees 1000gene-1000bp/estimated-genetrees 1000gene-true)
innerFolderNames48=(0.5X-1000-500 1X-25-500 1X-50-500 1X-100-500 1X-200-500 1X-500-500 1X-1000-500 2X-1000-500)
# outgroups=(11 O GAL STRCA)

# Rooting Methods MAD, MP, MV, OG, RD, NOCHANGE, RAND
methodNames=(  MP MV OG RAND )


# Inputs = Unrooted Gene trees, outputs = rooted gene trees to be used as stelar inputs

for folder in ${folders[@]}
do
    
	echo $folder
	speciesTree=../Dataset/$folder/true_tree_trimmed
	echo $speciesTree

	if [ $folder == "11-taxon" ];then
		innerFolderNames=("${innerFolderNames11[@]}")
		R=20
		outgroup='11'
	
	elif [ $folder == "15-taxon" ];then
		innerFolderNames=("${innerFolderNames15[@]}")
		R=10
		outgroup='O'
	
	elif [ $folder == "37-taxon" ];then
		innerFolderNames=("${innerFolderNames37[@]}")
		R=20
		outgroup='GAL'

	elif [ $folder == "48-taxon" ];then
		innerFolderNames=("${innerFolderNames48[@]}")
		R=20
		outgroup='STRCA'

	else
		echo ""
	fi
	echo  ${innerFolderNames[@]}

	
done



# Run stelar on rooted gene trees and get the average running time and RF distance
# Inputs = rooted gene trees, outputs = species tree

# Arrays to store sum and count for each method
declare -A sum_diffs
declare -A sum_rfs
declare -A count_diffs

# Global csv files for overall Summary
mkdir -p ../Dataset/Summary
rf_csv_main=../Dataset/Summary/average_rfs.csv

printf "%s,%s,%s,%s,%s\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv_main


# NEWW
rf_scores_file="rf_scores.txt"
echo -n "" > "$rf_scores_file" # Clear the file

# Function to calculate and store RF scores
calculate_and_store_rf() {
    folder=$1
    inner_folder=$2
    method=$3
    R=$4

    for ((j2=1; j2<=$R; j2++)); do
        stelar_output="../Dataset/$folder/$inner_folder/R$j2/stelar_outputs/stelar_output_$method.tre"
        
        if [ -f "$stelar_output" ]; then
            rf=$(python3 ../RF/getFpFn.py -e "$stelar_output" -t "$true_tree" | cut -d',' -f2 | tr -d '()')
            echo "$folder,$inner_folder,$method,$rf" >> "$rf_scores_file"
        fi
    done
}

# NEWW

bad=0
for folder in ${folders[@]}
do
	echo $folder
	
	# Local csv files for each taxa number
	mkdir -p ../Dataset/$folder/Summary
	rf_csv=../Dataset/$folder/Summary/average_rfs.csv
	printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv


	if [ $folder == "11-taxon" ];then
		innerFolderNames=("${innerFolderNames11[@]}")
		R=20
	
	elif [ $folder == "15-taxon" ];then
		innerFolderNames=("${innerFolderNames15[@]}")
		R=10
		echo "15 taxon"
		echo $R
	
	elif [ $folder == "37-taxon" ];then
		innerFolderNames=("${innerFolderNames37[@]}")
		R=20
	elif [ $folder == "48-taxon" ];then
		innerFolderNames=("${innerFolderNames48[@]}")
		R=20
	fi


	for inner_folder in ${innerFolderNames[@]}
	do
	 
		for method in  ${methodNames[@]}
		do
			# Initialize sum and count for each method
			sum_diffs[$method]=0
			count_diffs[$method]=0
			sum_rfs[$method]=0
		done
	    # temp_=0
	 
		for ((j=1;j<=$R;j++))
		do

			# echo "j at a: $j"
			mkdir -p ../Dataset/$folder/$inner_folder/R$j/stelar_outputs/
			# echo "j at b: $j"
			for method in  ${methodNames[@]}
			do
				# echo "j at c: $j"
				if [ ! -f "../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.tre" ]; then
					echo "all output files not present for $folder $inner_folder $j $method, current j: $j, current R: $R"
                    echo "skipping ********************************************************"
                    bad=1
                    break
				fi
				# echo "j at d: $j"
			    if [ $method == "RD" ]; then
						# echo "j at e: $j"
                        if [  "$inner_folder" == "100gene-true" ] || [ "$inner_folder" == "1000gene-true" ];then
                            echo "RD and skipping innerfolder $innerfolder"
							echo "j at f: $j"
							#
							j=1
                            # temp_=1
                            # break
                            continue
                        fi
				fi
				# temp_=0

				

				input=../Dataset/$folder/$inner_folder/R$j/stelar_inputs/stelar_input_$method.tre
				input_sanitized=../Dataset/$folder/$inner_folder/R$j/stelar_inputs/stelar_input_sanitized$method.tre
				output=../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.tre
				stripped_output=../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.stripped.tre
				true_tree=../Dataset/$folder/true_tree_trimmed
				# python3 ../RF/getFpFn.py -e ../Dataset/15-taxon/100gene-100bp/estimated-genetrees/R1/stelar_outputs/stelar_output_MP.tre -t ../Dataset/15-taxon/true_tree_trimmed

				echo $input $output $true_tree
				index_to_extract=2
				tuple=$(python3 ../RF/getFpFn.py -e $output -t $true_tree)
				RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
				echo $RFdistance

				# sum_diffs[$method]=$(echo "${sum_diffs[$method]} + $DIFF" | bc)
				sum_rfs[$method]=$(echo "${sum_rfs[$method]} + $RFdistance" | bc)
                ((count_diffs[$method]++))
				

				calculate_and_store_rf "$folder" "$inner_folder" "$method" "$R"
				# rm $input_sanitized.log
				# break

            done
			


		done

		

		# Print average DIFF for each method in local and global csv files

		for method in  ${methodNames[@]}
        do
            if [ ${count_diffs[$method]} -gt 0 ]; then

                average_rf=$(echo "${sum_rfs[$method]} / ${count_diffs[$method]}" | bc -l)
				printf "Average RF for %s: %f\n" "$method" "$average_rf"
                
                printf "stelar,%s,%s,%s,%.6f\n" "$folder" "$inner_folder" "$method" "$average_rf" >> $rf_csv
				printf "stelar,%s,%s,%s,%.6f\n" "$folder" "$inner_folder" "$method" "$average_rf" >> $rf_csv_main
              

            fi
        done
		
		printf "\n" >> $rf_csv
		printf "\n" >> $rf_csv_main

		# break
	done
	# break
done

python3 calculate_stats.py "$rf_scores_file" "rf_average.csv"

if [ $bad -eq 1 ];then
    echo "some files missing, incomplete run. Please ckeck stelar_outputs and rerun"
fi