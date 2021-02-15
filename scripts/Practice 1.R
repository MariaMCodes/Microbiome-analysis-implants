

#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

#Define a default theme for ggplot graphics
theme_set(theme_bw())


otu_table = read.csv("table-with-no-from_biom.csv", sep=",", row.names=1)
otu_table

tax_table = read.csv("data/taxonomy.csv", sep=",", row.names=1)
tax_table

metadata = read.csv("data/alcl_cc2.csv", sep=",", row.names=1)
metadata
#check what structure data is in
class(metadata)

#convert to matrix
OTU <- as.matrix(otu_table)
TAX <- as.matrix(tax_table)
#try if data is not in data.frame but because it already is, no need to do
#rownames(metadata) <- metadata[,1]  ## as data.frame
#samples <- as.matrix(metadata)      ## as matrix
#rownames(samples) <- metadata$Diagnosis


OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
samples = sample_data(metadata)
phy_tree = read_tree("data/rooted-tree.nwk")


#Make phyloseq object with matrix
phy <- phyloseq(OTU, TAX, samples, phy_tree)
phy

#to check row names are the same
rownames(metadata)
sample_names(phy)
#they are now the same

## PRACTICE 1 from website https://rforearth.wordpress.com/author/kerryleekuntz/
#plot bar plot
bar_plot <- plot_bar(phy, fill = "phylum")
bar_plot

#create phylogenetic tree
#random_tree = rtree(ntaxa(phy), rooted = TRUE, tip.label = taxa_names(phy))
