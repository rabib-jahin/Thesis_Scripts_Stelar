#!/usr/bin/perl

use lib "../bioPerl-1.5.2/lib/perl5";
use strict;
use warnings;
use Bio::TreeIO;



# Check for correct number of arguments
die "Usage: $0 <newick_tree_file>\n" unless @ARGV == 1;

# Read the Newick tree file
my $tree_file = $ARGV[0];

# Check if file is readable
unless (-r $tree_file) {
    die "Cannot read file $tree_file\n";
}

my $treeio = Bio::TreeIO->new(-file => $tree_file, -format => 'newick');

# Process the tree
if (my $tree = $treeio->next_tree) {
    my $root = $tree->get_root_node;
    if ($root && defined $root->id && $root->id ne '') {
        # Print only the root node ID
        print $root->id, "\n";
    } else {
        print "No named root found\n";  # or use 'die' to indicate an error
    }
} else {
    die "Error: Unable to parse the tree from $tree_file\n";
}

# ./get_tree_root.pl ../Dataset/11-taxon/estimated_Xgenes_strongILS/estimated_15genes_strongILS/R1/stelar_outputs/stelar_output_OG.tre
# ((Human, Chimpanzee), (Gorilla, (Orangutan, Gibbon)))Temp;

# ((Human, Chimpanzee), (Gorilla, (Orangutan, Gibbon)))Temp;
# (((Human, Chimpanzee), (Gorilla, (Orangutan, Gibbon))), Temp);