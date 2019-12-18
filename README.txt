#### READ ME #### 

######R code for: Parras-Moltó & Aguirre de Cárcer, 2019#######

##This document describes the use of FunCongr.R, as presented in “Assessment of phylo-functional coherence along the bacterial phylogeny and taxonomy” by Parras-Moltó and Aguirre de Cárcer (https://www.biorxiv.org/content/10.1101/795914v1).
##The procedure employed trees, taxonomic assignments, and 16S rRNA gene sequences from GTDB (https://gtdb.ecogenomic.org), as well as the full list of genome annotations from proGenomes (http://progenomes.embl.de/).

##Before the use of the script itself we prepared the input data; First we transformed genome annotations into a binary presence/absence table of gene content per genome (COG annotations were employed for consistency, representing >77% of the total annotations). 

##Then, for each pairwise comparison between genomes we produced 16S rRNA gene distances and Jaccard distances based on gene content (average of A-B and B-A), summarized in a four columns table (tab_info_pair.txt).

##Constant node labels were given to the input tree using makeNodeLabel from ape package in R (_name). Subsequently, the tree was pruned to contain only leaves with annotated genomes (_name_pruned).

##FunCongr.R reads the named_pruned tree, as well as the four columns table in the same directory, and prints out, fore each evaluated node, node name and hour (nodes passing the test) plus within-node average and standard deviation of Jaccard values (nodes failing the test, e.g. output_16s.log). The script does not evaluate nodes with fewer than 5 members or nodes presenting descendant nodes that already failed the test. Due to its nature, the script runs in a single processor, which, due to the size of the tree, translated into a lengthy (>400h) run.

##The results are later parsed with Paint_tree.R, which reads in the same directory the same tree and four columns table, as well as the metadata table containing the taxonomic affiliation of each leave (bac_metadata_r86.tsv) and a list of failing nodes identified with FunCongr.R (e.g. bad_nodes.txt). The scripts follows the rule that non-functionally coherent nodes are those that either failed the test or present a descendant node that failed the test. Then, for each last node were functional congruence was observed along each branching path it outputs a table with node ID, node label, number of leaves, and average plus standard deviation values of within-node 16S rRNA gene and Jaccard-based distances, and the node’s 80% consensus taxonomy) (nodes_info.txt), as well as a painted tree as a pdf image. 
