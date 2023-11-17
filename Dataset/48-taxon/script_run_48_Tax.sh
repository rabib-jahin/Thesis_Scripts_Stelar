#!/bin/bash

innerFolder_arr=('1X-1000-500/1X-50-500' '1X-1000-500/1X-100-500' '1X-1000-500/1X-200-500' '1X-1000-500/1X-500-500' 
	'0.5X-1000-500' '1X-1000-500' '2X-1000-500')

# here_folder_arr=('1X-25-500' '1X-50-500' '1X-100-500' '1X-200-500' '1X-500-500' '0.5X-1000-500' '1X-1000-500' '2X-1000-500')

## 25 genes is not shown ##

here_folder_arr=('1X-50-500' '1X-100-500' '1X-200-500' '1X-500-500' '0.5X-1000-500' '1X-1000-500' '2X-1000-500')

makeDirsAndCopyGTs()
{

	for (( iter = 0; iter < ${#innerFolder_arr[*]}; iter++ ))
	do
		folderName=${innerFolder_arr[$iter]}
		hereFolder=${here_folder_arr[$iter]}

		## make directories ##
		# mkdir "$hereFolder"

		for i in {1..20}; do
			
			## make directories ##
			# mkdir "$hereFolder/R$i"

			## copy gene-trees ##
			# cp "/home/mahim/Desktop/AVIAN_WHOLE/$folderName/R$i/all_gt.tre" "$hereFolder/R$i/all_gt.tre"

			## count lines ##
			# wc -l "$hereFolder/R$i/all_gt.tre"

			# diff "/home/mahim/Desktop/AVIAN_WHOLE/$folderName/R$i/weighted_quartets" "$hereFolder/R$i/weighted_quartets"

			# astral.5.7.3.heuristic.tre

			cp "/home/mahim/Desktop/AVIAN_WHOLE/$folderName/R$i/astral.5.7.3.heuristic.tre" "$hereFolder/R$i/astral-July26.5.7.3.tre"
		
		done
		echo
	done
}




makeDirsAndCopyGTs

