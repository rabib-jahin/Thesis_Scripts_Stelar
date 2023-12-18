#!/bin/bash
# set -ex
# java -jar STELAR.jar -i <input-gene-tree-name> -st <score-species-tree-name>


#!/bin/bash
# set -ex
# Local files

# Outer folders with taxa numbers: 11-taxon, 15-taxon, 37-taxon, 48-taxon
folders=( 11-taxon 15-taxon )

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

# # Path to your Perl script
# perl_script_path="./get_tree_root.pl"

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


# Global csv files for overall Summary
mkdir -p ../Dataset/Summary
# rf_csv_main=../Dataset/Summary/average_rfs.csv

# printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv_main


# NEWW
# rooting_matches_file="../Dataset/Summary/rooting_matches.csv"
# printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "method" "R" "rooting_match" > "$rooting_matches_file" # Clear the file
quartet_scores_file="../Dataset/Summary/quartet_scores.csv"
printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "method" "R" "quartet_score" > "$quartet_scores_file" # Clear the file


# echo -n "folder, inner-folder, method, R, rf-score" > "$rf_scores_file" # Clear the file

# NEWW

bad=0
for folder in ${folders[@]}
do
	echo $folder
	
	# Local csv files for each taxa number
	mkdir -p ../Dataset/$folder/Summary
	rf_csv=../Dataset/$folder/Summary/average_rfs.csv
	# printf "%s,%s,%s,%s,%s\n\n" "summary_method" "folder" "inner_folder" "rooting_method" "average_rf" > $rf_csv


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
				


				true_tree=../Dataset/$folder/true_tree_trimmed
				stelar_output="../Dataset/$folder/$inner_folder/R$j/stelar_outputs/stelar_output_$method.tre"
				# #get quartet score in quartet_score
				# # java -jar astral.5.7.8.jar -q test_data/simulated_14taxon.default.tre -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_scored.tre 2> test_data/simulated_scored.log
				# quartet_score=$(grep "Quartet score is:" "$astral_log" | awk '{print $4}')
                # echo "$folder,$inner_folder,$method,$j,$quartet_score" >> "$quartet_scores_file"
				# java -jar astral.5.7.8.jar -q test_data/simulated_14taxon.default.tre -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_scored.tre 2> test_data/simulated_scored.log

				astral_log="../Dataset/$folder/$inner_folder/R$j/astral_$method.log"
				java -jar astral.5.7.8.jar -q "$true_tree" -i "$stelar_output" -o /dev/null 2> "$astral_log"

				# Extract quartet score from log
				quartet_score=$(grep "Final quartet score is:" "$astral_log" | awk '{print $5}')
				echo "quartet score: $quartet_score"
				# Write quartet score to file
				echo "$folder,$inner_folder,$method,$j,$quartet_score" >> "$quartet_scores_file"
				# rm $astral_log
                
            done
			


		done

		
	done
	# break
done


if [ $bad -eq 1 ];then
    echo "some files missing, incomplete run. Please ckeck stelar_outputs and rerun"
fi