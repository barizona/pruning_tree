#!/usr/bin/env Rscript

# how to run from command line
# Rscript pruning_tree.R --original_tree gubbins.masked.snpsites.fasttree.nwk --alignment gubbins.masked.snpsites.aln

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Excluding too long branches ---------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxx
# Nextflow pipeline: after FastTree - before rooting_tree_R2_inputs_for_count.R
# Abbreviations:
#       IQR: the upper limit of the Inter Quartile Range of a vector: Q3 + 3 * (Q3 - Q1) = extreme outlier threshold
#       NT: the Number of Tips on the tree
# Aim: Excluding IQR branches having less than the 0.05 (or a given) proportion of the tips for unrooted trees.
#     Calculating the number of tips for each branch:
#            
#            Excluding a single branch: 
#                   Excluding a single top IQR branch (based on branch lengths) to avoid false midpoint rooting. 
#                   The branch can be excluded only if it has no more than the 0.05 (or a user defined proportion) of the tips.
#            Iterating till there are no more IQR tips based on the R2T distances:
#                  Midpoint rooting.
#                  Calculating R2T distances.
#                  Excluding the top 1 IQR tip based on the R2T distances.
#            Unrooting the pruned tree.
# Source files:
#    pruning_extreme_outlier_branches_and_tips.R
# Inputs: 
#    required:
#        original_tree: NEWICK tree file
#    optional:
#        alignment: sequence alignment FASTA file used for inferring the original tree; provide only if branch lengths should be multiplied by the length of the alignment
#        outgroup: name of the single outgroup to root the tree with
#        tipprop: proportion of tips on a single branch that can be excluded based on its IQR.
# Outputs: 
#    pruned.nwk: newick tree file
#    branch_length_distribution.png & pdf: histogram of branch lengths indicating the count of included and excluded tips
#    branch_length_nr_of_tips.png & pdf: scatter plot of branch lengths and number of tips with the IQR threshold and the number of tips threshold
#    original_vs_pruned_tree.png & pdf: tree plot indicating the pruned tip branches
#    tree_pruning.log: text file containing the summary comparison of the original tree and the pruned tree.
#xxxxxxxxxxxxxxxxxxxx

# Further name suggestions:
# Trimmed tree
# Pruned tree
# Reduced tree
# Filtered tree
# Cleaned tree
# Simplified tree
# Condensed tree
# Refined tree
# Revised tree
# Improved tree

#xxxxxxxxxxxxxxxxx
# Initialization ----------------------------------------------------------
#xxxxxxxxxxxxxxxxx

#xxxx
# * Call required packages -----
#xxxx
library(optparse) # named arguments
library(tidyverse)
library(magrittr) # %<>% pipe
library(treeio) # as_tibble
library(phytools) # midpoint.root
library(ggtree) # tree plotting
library(Biostrings) # readDNAStringSet at func_alignment_length function
library(caper) # clade.members.list at func_nr_of_tips function

#xxxx
# * Define arguments -----
#xxxx
option_list <- list(
  make_option(c("-t", "--original_tree"), 
              type = "character", default = NULL, 
              help = "input newick original_tree file from which long branches should be excluded", 
              metavar = "character"),
  make_option(c("-a", "--alignment"), 
              type = "character", default = NULL, 
              help = "sequence alignment used for inferring the original_tree", 
              metavar = "character"),
  make_option(c("-p", "--tipprop"), 
              type = "double", default = NULL, 
              help = "proportion of tips on a branch to drop if the branch is too long", 
              metavar = "character"))

opt_parser  <- OptionParser(option_list = option_list)

if (!interactive()) {
  opt  <- parse_args(opt_parser)
} else {
  opt  <- parse_args(opt_parser)
  opt$dir <- "../Main_STs_results/ST73/fasttree/"
  opt$original_tree <- paste0(opt$dir, "gubbins.masked.snpsites.fasttree.nwk")
  opt$alignment <- paste0(opt$dir, "gubbins.masked.snpsites.aln")
}


#xxxx
# * Delete the output files if exist -----
#xxxx
unlink(c(".RData",
         "pruned.nwk",
         "branch_length_distribution.png",
         "branch_length_distribution.png",
         "branch_length_nr_of_tips.png",
         "branch_length_nr_of_tips.pdf",
         "original_vs_pruned_tree.png",
         "original_vs_pruned_tree.pdf",
         "tree_pruning.log"))

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Functions -----------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxx
# Alignment length
#xxxxxxxxxxxxxxxxxx
# read alignments to get their lengths to multiply fasttree barnchlangths with it (-> SNPs / genome for Bactdate)
func_alignment_length <- function(fastafile) {
  alignment <- readDNAStringSet(fastafile, format = "fasta")
  # length of the first sequence
  a_length <- alignment[[1]]@length
  return(a_length)
}
