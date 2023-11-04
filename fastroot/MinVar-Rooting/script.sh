#!/bin/bash
# set -ex

timeFile="11-taxa-time.csv" # will keep the records of running time
RFRateFile="11-RFRateFile.csv" # will keep the records of RF rates


folders=( 11-taxon )
innerFolderNames=(estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS )
methodNames=( MP MV)
# Headers of the csv files
#echo "ilslevel, genetrees, basepair, methodName, modelConditionName,Robinson-Foulds distance" > $RFRateFile
#echo "name,time" > $timeFile

for folder in ${folders[@]}
do
	echo $folder
	speciesTree=$folder/true_tree_trimmed
	echo $speciesTree
	
	for inner_folder in ${innerFolderNames[@]}
	do
	  
	 
		for j in {1..20}
		do
			 #if [ $folder == "noscale.200g.true" ]; then
			  #gt=$folder/R$j/BS.1
  		 	#else
			#  gt=$folder/R$j/R$j.genetrees
			  #break
			# fi
			#echo $folder

			#sed 's/e-//g' $gt > $gt.temp #to replace e- with nothing
            #./strip_edge_support2.pl -i $gt.temp -o $gt.stripped
            #./reroot_tree.pl -t $gt.stripped -r GAL -o $gt.stripped.rooted
			
			gt=$inner_folder/R$j/all_gt.tre
			

			# if [[ -f "$speciesTree.stripped" ]];
			# then 
			# 	continue
			# fi
			for method in  ${methodNames[@]}

			do

			if [ $method == "MP" ]; then 
			      
					file=./$folder/$gt
					START=$(date +%s.%N)
					python3 utils.py $file temp.tre
					python3 FastRoot.py -m MP -i temp.tre  -o out.tree 
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
					# cp gt.best.of.10.tre  $speciesTree
					# rm -r gt*
					echo $DIFF
					rm temp.tre
					break
			
			else 		
			    break
				if [ $method == "SuperTriplets_v1.1" ]; then
					if [[ -f "$speciesTree" ]]; 
					then 
						rm -r $speciesTree
					fi
					START=$(date +%s.%N)
			 		java -jar $method.jar  $gt  $speciesTree
					END=$(date +%s.%N)
        	        DIFF=$(echo "$END - $START" | bc)
					

				else
					
					START=$(date +%s.%N)
					java -jar $method.jar -i $gt -o $speciesTree
					END=$(date +%s.%N)
					DIFF=$(echo "$END - $START" | bc)
				fi
			fi
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

