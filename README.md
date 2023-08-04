# IQR Tree Pruner

The aim of the *R* script is to exclude extremely long branches --
representing potentially misclassified taxa or genomic regions with
undetected recombinant parts -- from a phylogenetic tree.

## Running the code

```         
Rscript pruning_tree.R --original_tree tree.nwk
```

## Inputs and outputs

### Input file and argument

-   required:
    -   `--original_tree`: a NEWICK tree file
-   optional argument:
    -   `--tipprop`: The proportion of tips on a single branch that can
        be excluded based on the upper fence of the IQR of branch
        lengths.The default value is 0.05 (5%).

### Outputs files

-   **pruned.nwk**: the pruned NEWICK tree file

-   **branch_length_nr_of_tips.png** & **pdf**: a
    scatter plot of branch lengths and number of tips with the IQR upper
    fence and the proportion of tips gave by the `--tipprop` argument

-   **branch_length_distribution.png** & **pdf**: a
    histogram of branch lengths indicating the count of included and
    excluded tips

-   **original_vs_pruned_tree.png** & **pdf**: a
    tree plot indicating the pruned tip branches

-   **tree_pruning.log**: a text file containing the summary
    comparison of the original tree and the pruned tree.

## The IQR pruning algorithm

Abbreviations:

-   *IQR*: the upper fence of the *Inter Quartile Range* of a vector:
    Q3 + 3 \* (Q3 - Q1) = extreme outlier threshold

-   *R2T*: root-to-tip distance on a midpoint rooted tree

1<sup>st</sup> part: Pruning the unrooted tree based on the upper fence
of *IQR* branch lengths.

1.  Calculating the minimum number of tips for each branch (unrooted
    bifurcating tree, both directions are looked up, than the least tip
    number is chosen to represent the branch).

2.  Calculating the *IQR* for branch lengths.

3.  Excluding those extreme outlier branches containing less than 0.05
    (or a given) proportion of the tips.

2<sup>nd</sup> part: Pruning the midpoint rooted tree based on
root-to-tip distances.

Excluding tips based on the *R2T* *IQRs*, while iteratively midpoint
rooting and excluding the top gretaest extreme outlier.

4.  Midpoint rooting the tree (after pruning with method described in
    the 1<sup>st</sup> part).

5.  Calculating the *IQR* for root-to-tip distances.

6.  Excluding the most extreme outlier tip based on the *IQR* for
    root-to-tip distances.

7.  Repeating point 4. to 6. till there are no more extreme outlier
    *IQR* tip is found.

8.  Unrooting the pruned tree.
