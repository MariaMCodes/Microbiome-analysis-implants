
#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

#Define a default theme for ggplot graphics
theme_set(theme_bw())

#Import biom files
#rich_dense_biom = system.file("data/table-with-no-from_biom.csv", package="phyloseq")
#treefilename = system.file("data/rooted-tree.nwk", package = "phyloseq")
#refseqfilename = system.file("data/dna-sequences.fasta", package = "phyloseq")
#
#import_biom(rich_dense_biom, treefilename, refseqfilename, parseFunction = 
#            parse_taxonomy_greengenes)







otu_table = read.csv("data/table-with-no-from_biom.csv", sep=",", row.names=1)
otu_table

otu_matrix = read.csv("data/taxonomy.csv", sep=",", row.names=1)
otu_matrix = as.matrix(otu_matrix)

metadata = read.csv("data/alcl_cc.csv", sep=",", row.names=1)
metadata


OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(otu_matrix)
meta = sample_data(metadata)
phy_tree = read_tree("data/rooted-tree.nwk")


phyloseq_merged = phyloseq(OTU, TAX, meta, phy_tree)
phyloseq_merged
sample_names(phyloseq_merged)


plot_bar(phyloseq_merged, fill = "kingdom")
plot_bar

plot_bar(phyloseq_merged, fill = "")
plot_bar

plot_tree(phyloseq_merged, color = "genus", shape = "Diagnosis", size = "OTU_ID")
plot_tree

plot_richness(phyloseq_merged, x = "Diagnosis", colour = "Description")
