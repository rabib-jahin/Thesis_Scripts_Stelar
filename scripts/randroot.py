#!/usr/bin/env python3
import sys
import random
from ete3 import Tree

# Check command line arguments
if len(sys.argv) != 3:
    print("Usage: python3 randomly_root_trees.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        # Parse the tree
        tree = Tree(line.strip(), format=1)
        
        # Get all leaf nodes
        leaves = tree.get_leaves()
        
        # Select a random leaf node
        random_leaf = random.choice(leaves)
        
        # Set the tree's outgroup to the random leaf
        tree.set_outgroup(random_leaf)

        
        # Write the randomly rooted tree without branch lengths to the output file
        f_out.write(tree.write(format=9) + '\n')

print(f"Randomly rooted trees without branch lengths have been written to {output_file}")
