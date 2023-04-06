library(ggplot2)
library(ggtree)
library(ape)
library(tidyr)
library(dplyr)
library(pheatmap)


path <- "Desktop/influenza\ A/phylogenetic_tree_analysis/H1N1"
#path <- "Desktop/influenza\ A/phylogenetic_tree_analysis/H3N2"
setwd(path)
first_fs <- "14_15"
second_fs <- "15_16"

info <- read.csv("metadata.tsv", sep="\t", header=TRUE)
info$Year <- as.character(info$year)

file <- "tree_ha.nwk"

tree <- read.tree(file)
#tree <- root(tree, "A/Wellington/01/2004")
tree <- root(tree, "A/Iowa/01/2010")


distances <- cophenetic.phylo(tree)

only_first_fs <- info %>% filter(flu_season == first_fs)
only_second_fs <- info %>% filter(flu_season == second_fs)
rows_to_select <- match(only_first_fs$sequence_id, rownames(distances))
cols_to_select <- match(only_second_fs$sequence_id, colnames(distances))

selected_rows <- distances[rows_to_select, cols_to_select]

pheatmap(selected_rows)
p1 <- ggtree(tree) %<+% info +
  geom_tippoint(aes(color=Year)) +
  geom_treescale(x=0, y=45, fontsize=4, linesize=2, offset=2) +
  geom_tiplab(size=3) 
plot(p1)
