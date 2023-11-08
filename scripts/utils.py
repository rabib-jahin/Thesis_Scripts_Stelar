# import re
# import sys

# def add_branch_length(newick_str):

#     newick_with_lengths = re.sub(r'([^():,;]+)(?=[,);])', r'\1:1', newick_str)
#     newick_with_lengths = re.sub(r'(\))(?=\))', r'\1:1', newick_with_lengths)
#     newick_with_lengths = re.sub(r'(\))(?=,)', r'\1:1', newick_with_lengths)
#     return newick_with_lengths

# if len(sys.argv) != 3:
#     print("Usage: python script.py input_file output_file")
#     sys.exit(1)

# input_file = sys.argv[1]
# output_file = sys.argv[2]

# with open(input_file, 'r') as file:
#     trees = file.readlines()

# output_trees = []
# for tree in trees:
#     tree_with_lengths = add_branch_length(tree)
#     output_trees.append(tree_with_lengths)

# with open(output_file, 'w') as file:
#     file.writelines(output_trees)

import re
import sys

    
def add_branch_length(newick_str):
    if ":" in newick_str:
        return newick_str
    # Use a regular expression to add a length of 1 only if there is no existing length.
    newick_with_lengths = re.sub(r'([^():,;]+)(?=[,);])', r'\1:1', newick_str)
    newick_with_lengths = re.sub(r'(\))(?=\))', r'\1:1', newick_with_lengths)
    newick_with_lengths = re.sub(r'(\))(?=,)', r'\1:1', newick_with_lengths)
    return newick_with_lengths

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as file:
    trees = file.readlines()

output_trees = []
for tree in trees:
    # # Check if there's already a branch length specified and skip modification.
    # if re.search(r'(:\d+\.\d+|:\d+|\d+\.\d+|\d+)', tree):
	# #output_trees.append(tree)
	# pass
    # else:
    tree_with_lengths = add_branch_length(tree)
    output_trees.append(tree_with_lengths)

with open(output_file, 'w') as file:
    file.writelines(output_trees)
