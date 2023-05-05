#!/usr/bin/env Rscript

# how to run from command line
# Rscript pruning_tree.R --original_tree gubbins.masked.snpsites.fasttree.nwk --alignment gubbins.masked.snpsites.aln

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Excluding too long branches based on branch lengths and root to tip distances for unrooted trees ---------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxx
# Nextflow pipeline: after FastTree - before rooting_tree_R2_inputs_for_count.R
# Abbreviations:
#       IQR: the upper limit of the Inter Quartile Range of a vector: Q3 + 3 * (Q3 - Q1) = extreme outlier threshold
#       R2T: Root to tip distance of a midpoint rooted tree
# Aim: Excluding IQR branches having less than the 0.05 (or a given) proportion of the tips for unrooted trees.
#      Than excluding tips based on the R2T IQRs, while iteratively midpoint rooting and excluding the top longest outlier 
# Method:
#      Calculating the minimum number of tips for each branch (unrooted bifurcating tree, both directions are looked up).
#      Calculating IQR for branch lengths.
#      Iterating till there are no more IQR tips based on the R2T distances:
#          Midpoint rooting.
#          Calculating R2T distances.
#          Excluding the top 1 IQR tip based on the R2T distances.
#          Unrooting the pruned tree.
# Inputs: 
#    required:
#        original_tree: NEWICK tree file
#    optional:
#        alignment: sequence alignment FASTA file used for inferring the original tree; provide only if branch lengths should be multiplied by the length of the alignment
#        tipprop: proportion of tips on a single branch that can be excluded based on its IQR.
# Outputs: 
#    pruned.nwk: newick tree file
#    branch_length_nr_of_tips.png & pdf: scatter plot of branch lengths and number of tips with the IQR threshold and the number of tips threshold
#    branch_length_distribution.png & pdf: histogram of branch lengths indicating the count of included and excluded tips
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
         "branch_length_nr_of_tips.png",
         "branch_length_nr_of_tips.pdf",
         "branch_length_distribution.png",
         "branch_length_distribution.png",
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

#xxxxxxxxxxxxxx
# Change node numbers to node names
#xxxxxxxxxxxxxx
# nodenr_char: contains node numbers as character string
# tree_table: is a phylogenetic tree converted to tibble
func_nodenrs_to_nodenames <- function(nodenr_char, phy_tab){
  nodenr_char %>% 
    as.numeric() %>% 
    tibble(node = .) %>% 
    left_join(., phy_tab, by = "node") %>% 
    dplyr::select(label) %>% 
    pull()
}
# func_nodenrs_to_nodenames(nodenr_char = c(1, 2), phy_tab = tree_tab)

#xxxxxxxxxxxxxxxxxx
# Create the clade.members.list for tips and nodes of a rooted tree
#xxxxxxxxxxxxxxxxxx
func_nr_of_tips <- function(phy, phy_tab) {
  # list of tip names of all internal nodes (defined with nrs.)
  clade_members_list <- clade.members.list(phy = phy, 
                                           tip.labels = TRUE, 
                                           include.nodes = TRUE) %>% 
    # converting a list of lists to a simple list
    list_flatten()
  
  #xxxxxxxxxx
  # divide the list to 2 lists: tips & nodes
  #xxxxxxxxxx
  # Tips
  #xxxxx
  clade_members_list_tips <- clade_members_list[grep(x = names(clade_members_list), 
                                                     pattern = "_tips", 
                                                     value = TRUE)]
  # rename names: node nrs to labels
  names(clade_members_list_tips) %<>% 
    # delete the "_tips" part from the name
    gsub(x = ., pattern = "_tips", replacement = "") %>% 
    as.numeric() %>% 
    tibble(node = .) %>% 
    left_join(., phy_tab, by = "node") %>% 
    dplyr::select(label) %>% 
    pull()
  
  #xxxxx
  # Nodes (the list contains node nrs. and not node names)
  #xxxxx
  clade_members_list_nodes <- clade_members_list[grep(x = names(clade_members_list), 
                                                      pattern = "_tips", 
                                                      value = TRUE,
                                                      invert = TRUE)]
  rm(clade_members_list)
  
  # rename names: node nrs to labels
  names(clade_members_list_nodes) %<>%
    # delete the "_nodes" part from the name
    gsub(x = ., pattern = "_nodes", replacement = "") %>% 
    as.numeric() %>% 
    tibble(node = .) %>% 
    left_join(., phy_tab, by = "node") %>% 
    dplyr::select(label) %>% 
    pull()
  
  # change node nrs. to names in the list elements as well
  clade_members_list_nodes %<>% 
    lapply(., func_nodenrs_to_nodenames, phy_tab = phy_tab)
  
  # delete parts containing only a single node (the node itself as the list names)
  # example: 
  # $N263
  # [1] "N263"
  # -> delete elements which length is 1
  Nr_nodes <- sapply(clade_members_list_nodes, length)
  clade_members_list_nodes <- clade_members_list_nodes[names(Nr_nodes[Nr_nodes > 1])]
  rm(Nr_nodes)
  
  # delete node names that are the same as the list names
  # example:
  # $N261
  # [1] "N261" "N262" -> delete "N261" from the vector
  # -> delete the first element in each vector
  clade_members_list_nodes %<>% 
    lapply(., function(x) x[-1])
  
  # create a list of lists to return
  clade_members_list <- list(tip = clade_members_list_tips, 
                             node = clade_members_list_nodes)
  return(clade_members_list)
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

#xxxxxxxxxxxxxxxxxx
# Calculating the Root to tip distance of a rooted original_tree
#xxxxxxxxxxxxxxxxxx
func_R2T <- function(phy, phy_tab) {
  # delete "Root to tip distance" if it is presented in the tab
  phy_tab$`Root to tip distance` <- NULL
  # calculate new "Root to tip distance"
  phy_tab <- node.depth.edgelength(phy) %>% 
    as_tibble() %>%
    dplyr::rename("Root to tip distance" = "value") %>% 
    # add tip and node names
    bind_cols(label = c(phy$tip.label, phy$node.label), .) %>% 
    # add collection date
    left_join(phy_tab, ., by = "label")
  return(phy_tab)
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

# tree table containing is.tip column
tree_tab <- original_tree %>% 
  as_tibble() %>%
  # adding a column: is.tip
  mutate(is.tip = case_when(label %in% original_tree$tip.label ~ "Yes",
                            TRUE ~ "No"))

#xxxxxxx
# * Create clade_members_list containing descendants of all nodes ----
#xxxxxxx
# Adding Nr. of tips
clade_members_list <- func_nr_of_tips(phy = original_tree, 
                                      phy_tab = tree_tab)

# Nr of tips
Nr_tips <- sapply(clade_members_list$tip, length)

# Add a N1, N2 and `Nr. of tips` column to tree_tab
tree_tab %<>% 
  # add tip nrs.
  left_join(., enframe(Nr_tips, name = "label", value = "N1"), 
            by = "label") %>% 
  # add 1 as tip nr. when it is a tip
  mutate(N1 = case_when(is.na(N1) ~ 1,
                        TRUE ~ N1)) %>% 
  # Based on nr. of tips on the tree, calculate the nr. of tips to the other direction
  mutate(N2 = length(original_tree$tip.label) - N1) %>%
  # keep the minimum value
  rowwise() %>%
  mutate(`Nr. of tips` = min(N1, N2)) %>% 
  dplyr::select(-N1, -N2)

rm(Nr_tips)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Pruning based on branch lengths ---------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxx
# * Nr. of tips threshold ----
#xxxxxx
# Nr. of tips < tipprop or 0.05 of all tips
if(is.null(opt$tipprop)) {
  opt$tipprop <- 0.05
}

nr_tips_threshold <- length(original_tree$tip.label) * opt$tipprop

#xxxxxx
# * Branch lengths threshold ----
#xxxxxx
bl_threshold <- tree_tab %>% 
  func_IQR_upper_threshold(phy_tab = ., 
                           col_name = "branch.length")

#xxxxxx
# * Excluding branches ----
#xxxxxx

# excluding those branches where length > bl_threshold & Nr. of tips < nr_tips_threshold

# vector of tips and internal nodes to exclude
exc_bl <- tree_tab %>% 
  filter(branch.length > bl_threshold & 
         `Nr. of tips` < nr_tips_threshold) %>% 
  dplyr::select(label) %>% 
  pull()

# get the tips and further nodes based on exc_bl
exc_bl <- clade_members_list$node[exc_bl] %>% 
  unlist(use.names = FALSE) %>% 
  c(exc_bl, .) %>% 
  unique()

exc_bl <- clade_members_list$tip[exc_bl] %>% 
  unlist(use.names = FALSE) %>% 
  c(exc_bl, .) %>% 
  unique()

# drop these tips to create a shrunken tree
pruned_tree_bl <- drop.tip(phy = original_tree, tip = exc_bl)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Pruning based on root to tip distances of midpoint rooted tree --------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# apply midpoint rooting
pruned_tree_r2t <- pruned_tree_bl %>% 
  midpoint.root()

# create tibble from original_tree
pruned_tree_tab <- pruned_tree_r2t %>% 
  as_tibble() %>% 
  # adding a column: is.tip
  mutate(is.tip = case_when(label %in% pruned_tree_r2t$tip.label ~ "Yes",
                            TRUE ~ "No"))

# original R2T distances
pruned_tree_tab %<>% func_R2T(phy = pruned_tree_r2t, phy_tab = .)

# Excluding the tip having the greatest R2T distance longer than the IQR threshold
# rerooting, recalculating till there are no more outliers detected

tips2exclude <- pruned_tree_tab %>% 
  filter(is.tip == "Yes") %>% 
  arrange(desc(`Root to tip distance`)) %>% 
  filter(`Root to tip distance` > func_IQR_upper_threshold(phy_tab = ., 
                                                           col_name = "Root to tip distance")) %>% 
  dplyr::select(label) %>% 
  pull() 

while(length(tips2exclude) > 0) {
  pruned_tree_r2t  %<>% drop.tip(., tips2exclude[1]) # apply midpoint rooting if outgroup was not provided
  pruned_tree_r2t %<>% midpoint.root(.)
  # create tibble from original_tree
  pruned_tree_tab <- pruned_tree_r2t %>% 
    as_tibble() %>% 
    # adding a column: is.tip
    mutate(is.tip = case_when(label %in% pruned_tree_r2t$tip.label ~ "Yes",
                              TRUE ~ "No"))
  # R2T distance
  pruned_tree_tab %<>% func_R2T(phy = pruned_tree_r2t, phy_tab = .)
  # recalculate tips2exclude
  tips2exclude <- pruned_tree_tab %>% 
    filter(is.tip == "Yes") %>% 
    arrange(desc(`Root to tip distance`)) %>% 
    filter(`Root to tip distance` > func_IQR_upper_threshold(phy_tab = .,
                                                            col_name = "Root to tip distance")) %>% 
    dplyr::select(label) %>% 
    pull() 
}

# unroot tree
pruned_tree_r2t %<>%
  unroot()

# write pruned tree
pruned_tree_r2t %>% write.tree("pruned.nwk")

#xxxxxxxx
# * Create Excluded columns to tree tab for plotting ----
#xxxxxxxx

# vector of tips and internal nodes to exclude
exc_r2t <- c(setdiff(pruned_tree_bl$tip.label, pruned_tree_r2t$tip.label),
             setdiff(pruned_tree_bl$node.label, pruned_tree_r2t$node.label))

# get the tips and further nodes based on exc_r2t
exc_r2t <- clade_members_list$node[exc_r2t] %>% 
  unlist(use.names = FALSE) %>% 
  c(exc_r2t, .) %>% 
  unique()

exc_r2t <- clade_members_list$tip[exc_r2t] %>% 
  unlist(use.names = FALSE) %>% 
  c(exc_r2t, .) %>% 
  unique() %>% 
  sort()

# Add a new columns to tree_tab to define which tips & nodes were excluded
tree_tab %<>% 
  mutate(Excluded = case_when(label %in% c(exc_bl, exc_r2t) ~ "Yes",
                              TRUE ~ "No")) %>% 
  mutate(`Excluded tip` = case_when(is.tip == "Yes" &Excluded == "Yes"  ~ "Yes",
                                    TRUE ~ "No")) %>% 
  # as factor
  mutate(`Excluded tip` = factor(`Excluded tip`, levels = c("Yes", "No")))

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plotting ----------------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxx
# * Scatter plot of branch lengths and number of tips ----
#xxxxxx

# creating a tree_tab for plotting
tree_tab_plotting <- tree_tab %>% 
  # if branch.length is NA -> 0
  mutate(branch.length = case_when(is.na(branch.length) ~ 0,
                                   TRUE ~ branch.length))
# colour by density
tree_tab_plotting$Density <- func_density(tree_tab_plotting$branch.length, 
                                                   tree_tab_plotting$`Nr. of tips`, 
                                                   n = 100)
# scatterplot
p <- tree_tab_plotting %>% 
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

# saving the plot to png and pdf
ggsave(plot = p, filename = "branch_length_nr_of_tips.png",
       width = 6.27, height = 3.77)
ggsave(plot = p, filename = "branch_length_nr_of_tips.pdf",
       width = 6.27, height = 3.77)

#xxxxxx
# * Histogram of branch lengths ----
#xxxxxx

# plotting the distribution
p <- tree_tab %>% 
  filter(is.tip == "Yes") %>% 
  ggplot(aes(x = branch.length, fill = `Excluded tip`)) +
  geom_histogram() +
  xlab("Branch lengths") +
  scale_fill_manual(values = c("Yes" = "firebrick1", "No" = "black")) +
  theme_minimal()

# saving the plot to png and pdf
ggsave(plot = p, filename = "branch_length_distribution.png",
       width = 6.27, height = 3.77)
ggsave(plot = p, filename = "branch_length_distribution.pdf",
       width = 6.27, height = 3.77)

#xxxxxx
# * Tree plot ----
#xxxxxx
# it looses the branchlength when starting from tree_tab (??????????)
p <- original_tree %>%
  as_tibble() %>%
  mutate(Excluded = case_when(label %in% c(exc_bl, exc_r2t) ~ "Yes",
                              TRUE ~ "No")) %>% 
  as.treedata() %>% 
  ggtree(aes(color = Excluded), 
         layout = "equal_angle",
         size = 0.5) +
  scale_color_manual(values = c("Yes" = "firebrick1", "No" = "black")) +
  theme(legend.position = "bottom")

# saving the tree plot to png and pdf
ggsave(plot = p, filename = "original_vs_pruned_tree.png", 
       width = 8, height = 8)
ggsave(plot = p, filename = "original_vs_pruned_tree.pdf", 
       width = 8, height = 8)
