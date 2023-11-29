#!/bin/bash

# raxmlHPC -m GTRGAMMA -n best -s ../MSA/sim.11taxon/R1/bins.25.bin.1/bins.25.bin.1.xml -N 10 -p 12345
# raxmlHPC -m GTRGAMMA -n best -s ../MSA/15-taxon/1000gene-1000bp/estimated-genetrees/R1/1/1.fasta -N 10 -p 12345
# raxmlHPC -m GTRGAMMA -n best -s ../MSA/15-taxon/1000gene-100bp/estimated-genetrees/R10/1/1.fasta -N 10 -p 12345

# MSA/15-taxon/1000gene-100bp/estimated-genetrees/R10/1/1.fasta



# --------------------------------------------
# GT from Seq
# raxmlHPC -m GTRGAMMA -n best -s ../MSA/10-taxon/lower-ILS/estimated-genetrees/R1/1/truegene.fasta -N 10 -p 12345
# Remove the specific RAxML files for each run
# for i in {0..9}; do
#     rm -f "RAxML_log.best.RUN.$i"
#     rm -f "RAxML_parsimonyTree.best.RUN.$i"
#     rm -f "RAxML_result.best.RUN.$i"
# done

# # Remove the RAxML_info.best file
# rm -f "RAxML_info.best"

# echo "Temporary files have been deleted."

# --------------------------------------------






# --------------------------------------------

# MAD trial
# file='../Dataset/37-taxon/noscale.25g.500b/R1/all_gt.tre'
# file='../Dataset/11-taxon/estimated_Xgenes_strongILS/estimated_5genes_strongILS/R1/all_gt.tre'

# # python3 ../mad/mad.py $file -n

# python3 ../scripts/utils.py $file temp.tre
# python3 ../mad/mad.py temp.tre -n

# --------------------------------------------
# stelar trial
# file=../Dataset/37-taxon/noscale.25g.500b/R1/all_gt.tre
# java -jar ../scripts/stelar.jar -i $file -o myout.tre	

# --------------------------------------------
# RF score trial


# --------------------------------------------