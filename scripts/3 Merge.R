### SEQUENCING ANALYSIS USING PHYLOSEQ
##Merge data

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")
theme_set(theme_bw())

otu_table = read.csv("data/table-silva_biom.csv", sep=",", row.names = 1)
otu_table

tax_table = read.csv("data/taxonomy.csv", sep=",", row.names=1)   #so R will automatically convert all the “NULL” entries in the dataset into NA
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
#otu_table()   OTU Table:         [ 1341 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1341 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1341 tips and 1323 internal nodes ]

phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]
#lost 4 taxa from the OTU table, meaning out amplification seems to have worked well to only target bacteria

#created new dataset name to keep it separate from section 2
BIm <- phy_sub
BIm 

#1)add a new sample variable to separate human samples from controls
humantypes = c("ALCL", "CC", "Contralateral", "Implant exchange", "Ruptured")
sample_data(BIm)$Human <- get_variable(BIm, "Diagnosis") %in% humantypes

BIm
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]

sample_variables(BIm)
#[1] "Diagnosis"        "ProcessingMethod" "Human"  

#1) Merge
mergedBI = merge_samples(BIm, "Diagnosis")
SD = merge_samples(sample_data(BIm), "Diagnosis")
print(SD[, "Diagnosis"])
print(mergedBI)
sample_names(BIm)
sample_names(mergedBI)
identical(SD, sample_data(mergedBI))


sample_data(mergedBI)$Diagnosis = sample_names(mergedBI)
sample_data(mergedBI)$Human = sample_names(mergedBI) %in% humantypes
plot_richness(mergedBI, "Human", "Diagnosis", title = "merged")


plot_richness(BIm, x = "Diagnosis", measures = c("Chao1", "Shannon"))
