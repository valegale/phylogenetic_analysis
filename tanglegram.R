library(ape)
library(phytools)
library(dendextend)
library(viridis)
library(dplyr)
library(phylogram)

path <- "Desktop/influenza\ A/phylogenetic_tree_analysis"
setwd(path)

tree1 <- read.tree(file = "tree_ha.nwk")
tree1 <- root(tree1, "A/Wellington/01/2004")
#is.ultrametric(tree1_ultrametric)


dnd1 <- as.dendrogram(tree1)

tree2 <- read.tree(file = "tree_na.nwk")
tree2 <- root(tree2, "A/Wellington/01/2004")

dnd2 <- as.dendrogram(tree2)

## rearrange in ladderized fashion
dnd1 <- ladder(dnd1)
dnd2 <- ladder(dnd2)
## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 5)

dndlist %>% ladderize %>% 
  # untangle(method = "step1side", k_seq = 3:20) %>%
  set("rank_branches") %>%
  tanglegram(common_subtrees_color_branches = TRUE, margin_inner = 5)
