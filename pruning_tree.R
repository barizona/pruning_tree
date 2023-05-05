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
  opt$dir <- "../../Main_STs_results/ST73/fasttree/"
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

#xxxxxxxxxxxxxxxxx
# Convert phy
#xxxxxxxxxxxxxxxx
# read trees calculated by Nextflow
func_convert_phy <- function(phy) {
  
  # transforms all multichotomies into a series of dichotomies
  phy %<>% multi2di() %>% 
    # unroot the tree
    unroot()
  
  # multiply the branch length by the alignment length if alignment file was provided
  if(!is.null(opt$alignment)) {
    phy$edge.length <- phy$edge.length * func_alignment_length(opt$alignment)
  }
  
  # name nodes if there are no node labels provided
  if(anyNA(as.numeric(phy$node.label))) {
    # Naming nodes
    phy$node.label <- paste0("N", 1:length(phy$node.label))
  }
  
  return(phy)
}

#xxxxxxxxxxxxxxxxxx
# Create tree table, add is.tip column and the minimum number of descendant tips
#xxxxxxxxxxxxxxxxxx

func_phy_tab <- function(phy) {
  # create tibble from original_tree
  phy_tab <- phy %>% 
    treeio::as_tibble() %>% 
    # adding a column: is.tip
    mutate(is.tip = case_when(label %in% original_tree$tip.label ~ "Yes",
                              TRUE ~ "No"))
  
  # list of tip names of all internal nodes (defined with nrs.)
  clade_members_list <- clade.members.list(phy = phy, 
                                           tip.labels = TRUE, 
                                           include.nodes = TRUE) %>% 
    # converting a list of lists to a simple list
    list_flatten()
  
  # calculate the Nr. of tips for each branch to one direction
  phy_tab$N1 <- sapply(phy_tab$node, function(x) {
    index <- which(phy_tab$node == x)
    # if node is a tip, the number of descendant tips is 1.
    if (phy_tab$is.tip[index] == "Yes") {
      return(1) 
      # else the number of descendant tips is counted from clade_members_list
    } else {
      index2 <- which(names(clade_members_list) == paste0(x, "_tips"))
      if (length(index) == 0) {
        msg <- paste0(
          "Failed to count descendant tips for node: ",
          x
        )
        stop(msg)
      }
      length(clade_members_list[[index2]])
    }
  })
  
  # add this 
  # Based on nr. of tips on ther tree, calculate the nr. of tips to the other direction
  phy_tab %<>% 
    mutate(N2 = length(original_tree$tip.label) - N1) %>% 
    # keep the minimum value 
    rowwise() %>%
    mutate(`Nr. of tips` = min(N1, N2))
  
  # delete columns
  phy_tab$N1 <- NULL
  phy_tab$N2 <- NULL
  
  return(phy_tab)
}

#xxxxxxxxxxxxxxxxxx
# IQR outlier upper threshold
#xxxxxxxxxxxxxxxxxx

func_IQR_upper_threshold <- function(phy_tab, col_name) {
  # Quartiles of R2T distances
  qtr_vec <- phy_tab %>% 
    filter(!is.na(!!sym(col_name)) &
             !!sym(col_name) > 0) %>% 
    dplyr::select(!!sym(col_name)) %>% 
    pull() %>% 
    quantile()
  # Inter Quartile Range of R2T distances: upper limit Q3 + 3 * (Q3 - Q1) - threshold
  R2T_threshold <- qtr_vec["75%"] + 3 * (qtr_vec["75%"] - qtr_vec["25%"])
  names(R2T_threshold) <- NULL
  return(R2T_threshold)
}

#xxxxxxxxxxxx
# Density function for colouring a scatter plot by density
#xxxxxxxxxxxx
func_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Inputs ---------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# tree
original_tree <- read.tree(opt$original_tree) %>% 
  func_convert_phy()

# tree table containing Nr. of tips
original_tree_tab <- func_phy_tab(original_tree)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Thresholds --------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxx
# * Tips ----
#xxxxxx
# Nr. of tips < tipprop or 0.05 of all tips
if(is.null(opt$tipprop)) {
  opt$tipprop <- 0.05
}

nr_tips_threshold <- length(original_tree$tip.label) * opt$tipprop

#xxxxxx
# * Branch lengths ----
#xxxxxx
bl_threshold <- original_tree_tab %>% 
  func_IQR_upper_threshold(phy_tab = ., 
                           col_name = "branch.length")

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Pruning -----------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Excluding all tips having greater R2T distances than the IQR threshold, except the outgroup tip
pruned_tree <- original_tree_tab %>% 
  filter(is.tip == "Yes" & 
           `branch.length` > bl_threshold &
           label != opt$outgroup) %>% 
  select(label) %>% 
  pull() %>% 
  drop.tip(original_tree, .)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxx
# * Scatter plot of branch lengths and number of tips ----
#xxxxxx

# creating a tree_tab for plotting
original_tree_tab_plotting <- original_tree_tab %>% 
  # if branch.length is NA -> 0
  mutate(branch.length = case_when(is.na(branch.length) ~ 0,
                                   TRUE ~ branch.length))
# colour by density
original_tree_tab_plotting$Density <- func_density(original_tree_tab_plotting$branch.length, 
                                                  original_tree_tab_plotting$`Nr. of tips`, 
                                                  n = 100)
# scatterplot
original_tree_tab_plotting %>% 
  ggplot(aes(x = branch.length, y = `Nr. of tips`,
             color = Density)) +
  geom_point() +
  viridis::scale_color_viridis(option = "turbo") +
  xlab("Branch length") +
  # add branchlength threshold
  geom_vline(aes(xintercept = bl_threshold, 
                 linetype = "Branch length outliers\nto the right"), 
             color = "darkred") +
  # add nr. of tips threshold
  geom_hline(aes(yintercept = nr_tips_threshold, 
                 linetype = paste0("Nodes with less than ", opt$tipprop*100, "%\nof the tips are below")),
             color = "darkgreen") +
  # legend for lines
  scale_linetype_manual(name = "Outlier detection", values = c(2, 2), 
                        guide = guide_legend(override.aes = 
                                               list(color = c("darkred", 
                                                              "darkgreen")))) +
  theme_minimal()

