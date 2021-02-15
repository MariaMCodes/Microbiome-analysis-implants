



#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

#Define a default theme for ggplot graphics
theme_set(theme_bw())

#import data. Note: items within data eg sample names need to match. Also had to convert data from
#biom files in csv files in conda as I had issues reading the biom files into R
otu_table = read.csv("table-with-no-from_biom.csv", sep=",", row.names=1)
otu_table

tax_table = read.csv("data/taxonomy2.csv", sep=",", row.names=1)
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
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1270 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1270 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1270 tips and 1253 internal nodes ]

#double check row names are the same - needed to add 'X' in front of sample ID and change the hiphen to a dot as biom files had this
rownames(metadata)
sample_names(phy)
#they are now the same


## PRACTICE 3 from website https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html
#filter out everything but the kingdom of Bacteria, since it is not our intent to amplify anything else
phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1266 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1266 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1266 tips and 1249 internal nodes ]
#lost 4 taxa from the OTU table, meaning out amplification seems to have worked well to only target bacteria

#rename data to track easily
P1 <- phy_sub
P1 

tree_phylum <- plot_tree(P1, color = "phylum", shape = "Diagnosis", size = "Abundance")
tree_phylum

tree_kingdom <- plot_tree(P1, color = "kingdom", shape = "Diagnosis", size = "Abundance")
tree_kingdom

bar_phylum <- plot_bar(P1, fill = "phylum")
bar_phylum

bar_class <- plot_bar(P1, fill = "class")
bar_class

bar_genus <- plot_bar(P1, fill = "genus")
bar_genus


#Preprocessing
#filtering
P1r  = transform_sample_counts(P1, function(x) x / sum(x) )
P1fr = filter_taxa(P1r, function(x) var(x) > 1e-5, TRUE)

#subsetting
P1.fir = subset_taxa(P1, phylum=="Firmicutes")
P1.fir = prune_samples(sample_sums(P1.fir)>=20, P1.fir)

#merging
P1.fir.merged = merge_taxa(P1.fir, taxa_names(P1.fir)[1:5])

#agglomerate the "Bacteroidetes-only" dataset at the taxonomic rank of family, and created an annotated tree of the result
p1sfb = subset_taxa(P1fr, phylum == "Bacteroidetes")
p1sfbg = tax_glom(p1sfb, "family")
plot_tree(p1sfbg, color="Diagnosis", shape="class", size="Abundance")
#agglomerate the "Proteobacteria-only" (pseudomonas) dataset at the taxonomic rank of family, and created an annotated tree of the result
p1sfp = subset_taxa(P1fr, phylum == "Proteobacteria")
p1sfpg = tax_glom(p1sfp, "family")
plot_tree(p1sfpg, color="Diagnosis", shape="class", size="Abundance")
#agglomerate the "Firmicutes-only" (staphylococcus) dataset at the taxonomic rank of family, and created an annotated tree of the result
p1sff = subset_taxa(P1fr, phylum == "Firmicutes")
p1sffg = tax_glom(p1sff, "family")
plot_tree(p1sffg, color="Diagnosis", shape="class", size="Abundance")

#transform abundance values to fractional abundance
transform_sample_counts(P1.fir, function(OTU) OTU/sum(OTU) )
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 338 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 338 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 338 tips and 327 internal nodes ]


#remove taxa occurring less than 3 times in at least 20% of samples. This protects against an OTU with small mean and trivially large CV
P2 = filter_taxa(P1, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

#define a human versus non-human categorical variable, and add this new variable to sample data
sample_data(P2)$human = factor(get_variable(P2, "Diagnosis") %in% c("ALCL", 
    "Contralateral", "CC", "Ruptured", "Implant exchange", "Negative control", "Positive control"))

#standardize abundances to the median sequencing depth
total = median(sample_sums(P2))
standf = function(x, t=total) round(t * (x / sum(x)))
p1s = transform_sample_counts(P2, standf)

#filter the taxa using a cutoff of 3.0 for the coefficient of variation
p1sf = filter_taxa(p1s, function(x) sd(x)/mean(x) > 3.0, TRUE)

#subset the data to Bacteroidetes, used in some plots
p1sfb = subset_taxa(p1sf, phylum=="Bacteroidetes")
#subset the data to Proteobacteria, used in some plots
p1sfp = subset_taxa(p1sf, phylum = "Proteobacteria")
#subset the data to Firmicutes, used in some plots
p1sff = subset_taxa(p1sf, phylum = "Firmicutes")

#graphic summary of Bacteroidetes - summarise this slice of the data with some graphics
title = "Bacteroidetes-only"
plot_bar(p1sfb, "Diagnosis", "Abundance", title=title)
#add Family level
plot_bar(p1sfb, "Diagnosis", "Abundance", "family", title=title)
#facetted into diagnosis type
plot_bar(p1sfb, "family", "Abundance", "family", 
         title=title, facet_grid="Diagnosis~.")

#graphic summary of Proteobacteria - summarise this slice of the data with some graphics
title = "Proteobacteria-only"
plot_bar(p1sfp, "Diagnosis", "Abundance", title=title)
#add Family level
plot_bar(p1sfp, "Diagnosis", "Abundance", "family", title=title)
#facetted into diagnosis type
plot_bar(p1sfp, "family", "Abundance", "family", 
         title=title, facet_grid="Diagnosis~.")

#graphic summary of Firmicutes - summarise this slice of the data with some graphics
title = "Firmicutes-only"
plot_bar(p1sff, "Diagnosis", "Abundance", title=title)
#add Family level
plot_bar(p1sff, "Diagnosis", "Abundance", "family", title=title)
#facetted into diagnosis type
plot_bar(p1sff, "family", "Abundance", "family", 
         title=title, facet_grid="Diagnosis~.")
