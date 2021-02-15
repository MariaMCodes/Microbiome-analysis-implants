

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


## PRACTICE 2 from website https://rpubs.com/faysmith/metabarcoding
#1) filter out everything but the kingdom of Bacteria, since it is not our intent to amplify anything else
phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub
#> phy
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1270 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1270 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1270 tips and 1253 internal nodes ]

#> phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1266 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1266 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1266 tips and 1249 internal nodes ]

#from output, only lost 4 taxa from the OTU table, meaning our amplification seems to have worked well to only target bacteria

#2) Normalise the OTUs to account for sampling/sequence differences
#There is some randomization at work here, so better set some seed (so the randomized numbers are reproducable) rarefy_even_depth is a function that will randomly 
#sample OTUs from each group until it reaches the same number of OTUs that is in the smallest sample size. There are some that claim that rarefaction should not be 
#used on this kind of data, but everyone does it anyway. Some fancier methods of normalizing the data exist, but they are more complicated and not as well published on.
set.seed(100)
#randomly sample to the lowest denomination
phy_sub <- rarefy_even_depth(phy_sub, rngseed = TRUE)
phy_sub
#> phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 990 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 990 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 990 tips and 974 internal nodes ]
#276OTUs were removed because they are no longer present in any sample after random subsampling 
#went from 1266OTUs in our dataset to 990 after normalizing for sample size



#plot the phylum composition of patients
library(ggplot2)
library(viridis)
library(hrbrthemes)

phy_phylum <- phy_sub %>%
  tax_glom(taxrank = "phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt () %>%
  filter(Abundance > 0.02) %>%
  arrange(phylum)
phy_phylum

ggplot(phy_phylum, aes(x = Sample, y = Abundance, fill = phylum)) +
  facet_grid("Diagnosis~.") +
  geom_bar(position="stack", stat = "identity") +
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
             scale_x_discrete("Sample") +
             guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
             ylab("Relative Abundance (Phyla > 2%) \n") +
             ggtitle("Phylum Composition of patients \n Bacterial communitites by Diagnosis") 
             






