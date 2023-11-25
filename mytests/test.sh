#!/bin/bash
# set -ex

file=../Dataset/11-taxon/estimated_Xgenes_strongILS/estimated_5genes_strongILS/R1/all_gt.tre

# python3 ../scripts/utils.py $file temp.tre
# python3 ../fastroot/MinVar-Rooting/FastRoot.py -m MV -i temp.tre  -o mv_rooted.tre
# python3 ../fastroot/MinVar-Rooting/FastRoot.py -m MV -i $file  -o mv_rooted.tre

# perl reroot_tree.pl -t mv_rooted.tre -r '11' -o ogrooted.tre
#!/bin/bash








perl reroot_tree.pl -t $file -r '11' -o temp.tre

# Usage: ./script.sh outgroup input_file output_file

# Assign command line arguments to variables
# outgroup="$1"
# input="$2"
# output="$3"

outgroup='11'
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
sed -e "s/$outgroup//g" -e "s/:0/,$outgroup/g" "$input" > "$output"

echo "Processing complete. Output written to $output"

java -jar ../scripts/stelar.jar -i ../mytests/stelar_input.tre -o ../mytests/stelar_output.tre



