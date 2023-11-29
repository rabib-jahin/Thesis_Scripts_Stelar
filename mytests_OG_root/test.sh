#!/bin/bash
# set -ex

file=../Dataset/11-taxon/estimated_Xgenes_strongILS/estimated_5genes_strongILS/R1/all_gt.tre
file1=../Dataset/37-taxon/noscale.25g.500b/R1/all_gt.tre
python3 ../scripts/utils.py $file temp.tre
# python3 ../fastroot/MinVar-Rooting/FastRoot.py -m MV -i temp.tre  -o mv_rooted.tre
# python3 ../fastroot/MinVar-Rooting/FastRoot.py -m MV -i $file  -o mv_rooted.tre

# perl reroot_tree.pl -t mv_rooted.tre -r '11' -o ogrooted.tre
#!/bin/bash


# outgroup='GAL'
# outgroup='11'





# perl reroot_tree.pl -t mv_rooted.tre -r $outgroup -o temp.tre
perl reroot_tree.pl -t $file -r $outgroup -o temp.tre


# Assign command line arguments to variables
# outgroup="$1"
# input="$2"
# output="$3"


input='temp.tre'
output="stelar_input.tre"

# Check if input file exists
if [ ! -f "$input" ]; then
    echo "Input file does not exist: $input"
    exit 1
fi

# Process the file
# -e: Execute script
# -i: Edit files in place (makes backup if extension supplied)
# sed -e "s/$outgroup//g" -e "s/:0)/,$outgroup)/g" "$input" > "$output"
sed -E "s/:[0-9]+(\.[0-9]+)?\)$outgroup;/,$outgroup);/g" "$input" > "$output"

echo "Processing complete. Output written to $output"

java -jar ../scripts/stelar.jar -i ../mytests/stelar_input.tre -o ../mytests/stelar_output.tre
# java -jar ../scripts/stelar.jar -i ../mytests/mv_rooted.tre -o ../mytests/stelar_output.tre



# ..\MSA\10-taxon\lower-ILS\estimated-genetrees\R1\1\truegene.fasta

# raxmlHPC-8.0.19-SSE3 -m GTRGAMMA -n best -s ..\MSA\10-taxon\lower-ILS\estimated-genetrees\R1\1\truegene.fasta -N 10 -p 12345
raxmlHPC -m GTRGAMMA -n best -s ../MSA/10-taxon/lower-ILS/estimated-genetrees/R1/1/truegene.fasta -N 10 -p 12345