

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

theme_set(theme_bw())

otu_table = read.csv("data/table-with-no-from_biom.csv", sep=",", row.names = 1)
otu_table

tax_table = read.csv("data/taxonomy.csv", sep=",", row.names=1, na = "")   #so R will automatically convert all the “NULL” entries in the dataset into NA
tax_table

metadata = read.csv("data/alcl_cc2.csv", sep=",", row.names=1)
metadata

OTU <- as.matrix(otu_table)
TAX <- as.matrix(tax_table)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
samples = sample_data(metadata)
phy_tree = read_tree("data/rooted-tree.nwk")

phy <- phyloseq(OTU, TAX, samples, phy_tree)
phy
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1270 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1270 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1270 tips and 1253 internal nodes ]

phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub