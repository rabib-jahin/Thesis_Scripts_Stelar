#!/usr/bin/env python3
from ete3 import Tree
import sys

def compare_roots(tree1_path, tree2_path):
    try:
        # tree1 = Tree(tree1_path)
        # tree2 = Tree(tree2_path)
        tree1 = Tree('(1,(2,3))')
        tree2 = Tree('(1,(2,3))')
        if tree1.get_tree_root().name == tree2.get_tree_root().name:
            return 1
        else:
            return 0
    except Exception as e:
        print(f"Error comparing trees: {e}", file=sys.stderr)
        return 0
    
compare_roots(sys.argv[1], sys.argv[2])