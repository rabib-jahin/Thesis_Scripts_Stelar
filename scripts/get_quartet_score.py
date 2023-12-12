#!/usr/bin/env python3
from ete3 import Tree
import itertools
import sys

def compare_quartet_trees(tree1, tree2, quartet):
    # Prune trees to contain only the quartet taxa
    pruned_tree1 = tree1.copy()
    pruned_tree2 = tree2.copy()
    
    pruned_tree1.prune(quartet, preserve_branch_length=True)
    pruned_tree2.prune(quartet, preserve_branch_length=True)

    # Compare topology
    return pruned_tree1.compare(pruned_tree2, unrooted=True)['norm_rf'] == 0

def get_quartet_score(tree1, tree2):
    taxa = tree1.get_leaf_names()
    matching_count = 0

    for quartet in itertools.combinations(taxa, 4):
        if compare_quartet_trees(tree1, tree2, quartet):
            matching_count += 1

    return matching_count

def main():
    if len(sys.argv) != 3:
        print("Usage: python quartet_score.py <tree1.newick> <tree2.newick>", file=sys.stderr)
        sys.exit(1)

    tree1_file, tree2_file = sys.argv[1], sys.argv[2]

    # Load trees
    tree1 = Tree(tree1_file)
    tree2 = Tree(tree2_file)

    # Calculate matching quartet count
    matching_quartets = get_quartet_score(tree1, tree2)
    print(matching_quartets)

if __name__ == "__main__":
    main()
