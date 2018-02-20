#!/usr/bin/env Rscript

# plot significant hits from treeWAS

# dependencies
library(ggtree)
library(dplyr)
library(phangorn)

# input
data_file <- "/home/sb2145/Desktop/BIGSdb/TreeWAS_Pangenome/C_coli_update/treeWAS/results/chicken.plot_data"
metadata_file <- "/home/sb2145/Desktop/BIGSdb/TreeWAS_Pangenome/C_coli_update/C_coli.representative_metadata.tsv"
tree_file <- "/home/sb2145/Desktop/BIGSdb/TreeWAS_Pangenome/C_coli_update/binary_presence_absence.nwk"
output_file <- "/home/sb2145/Desktop/BIGSdb/TreeWAS_Pangenome/C_coli_update/treeWAS/results/chicken_plot.pdf"

# read tree
tree <- midpoint(read.tree(tree_file))

# read plot data 
pd <- read.delim(data_file, sep = "\t")
pd_in <- pd %>% select(-id)
rownames(pd_in) <- pd$id

# read metadata
meta <- read.delim(metadata_file, sep = "\t")
meta_plot <- data.frame(id = meta[,1], Host = (meta$Host == "chicken")) 

# prepare tree plot
t <- ggtree(tree) %<+% meta_plot + 
  geom_tiplab(aes(colour=Host),size=1, align=T, linesize = 0.1) + 
  scale_color_manual(values=c("blue", "red"))  

# plot heatmap
gheatmap(t, pd_in, offset = 0.01, width=0.5, font.size=2, colnames_angle=-45, hjust=0) +
  theme(legend.position="bottom")
?gheatmap
