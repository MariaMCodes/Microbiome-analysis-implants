### SEQUENCING ANALYSIS USING PHYLOSEQ
## SECTION 1: Importing Data and Functions for Accessing Data

#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

#Define a default theme for ggplot graphics
theme_set(theme_bw())

#import data. Note: items within data eg sample names need to match. Also had to convert data from
#biom files in csv files in conda as I had issues reading the biom files into R
otu_table = read.csv("data/table-silva_biom.csv", sep=",", row.names = 1)
otu_table

tax_table = read.csv("data/taxonomy.csv", sep=",", row.names=1, na = "")   #so R will automatically convert all the “NULL” entries in the dataset into NA
tax_table

metadata = read.csv("data/alcl_cc2.csv", sep=",", row.names=1)
metadata
#check what structure data is in
class(metadata)

#try importing ref seq file but can't merge into phyloseq
rs_file = read.csv("data/dna-sequences.fasta")
rs_file

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
#otu_table()   OTU Table:         [ 1341 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1341 taxa by 7 taxonomic ranks ]          # only 6 taxanomic ranks but have now inclueded a Gram pos or Gram neg column
#phy_tree()    Phylogenetic Tree: [ 1341 tips and 1323 internal nodes ]

#double check row names are the same - needed to add 'X' in front of sample ID and change the hiphen to a dot as biom files had this
rownames(metadata)
sample_names(phy)
#they are now the same


## PRACTICE 3 from website http://joey711.github.io/phyloseq/import-data.html
# 1) Filter out everything but the kingdom of Bacteria, since it is not our intent to amplify anything else
phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 2 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]
#lost 74 taxa from the OTU table, meaning out amplification seems to have worked well to only target bacteria

#rename data to track easily (BI = breast implants)
BI <- phy_sub
BI 

# Accessor functions 
ntaxa(BI)
#[1] 1267

nsamples(BI)
#[1] 79

sample_names(BI) [1:10]
#[1] "X1612LC" "X1612RC" "X1618LC" "X1618RC" "X1620LC" "X1621LC" "X1621LI"
#[8] "X1621RC" "X1626LC" "X1626RC"

rank_names(BI)
#[1] "kingdom" "phylum"  "class"   "order"   "family"  "genus"  

sample_variables(BI)
#[1] "Diagnosis"        "ProcessingMethod"

otu_table(BI) [1:5, 1:5]                   # [rows, column]
#OTU Table:          [5 taxa and 5 samples]
#                    taxa are rows
#                                 X1612LC X1612RC X1618LC X1618RC X1620LC
#d0a102d533aa6170c2505bbef5ec7101       0       0       0       0       0
#7aaa3ecae849160810cc9fedc59df606       0       0       0       0       0
#ef9faac34d710ef9840e9a1afa1ee09f       0       0       0       0       0
#6f0967215a05e1eeb02ad466bea710bf       0       0       0       0       0
#b141f55e90b2532c91dabda3907f973e       0       0       0       0       0

tax_table(BI) [1:5, 1:5]
#Taxonomy Table:     [5 taxa by 5 taxonomic ranks]:
#                                 kingdom    phylum             class             
#d0a102d533aa6170c2505bbef5ec7101 "Bacteria" "Cyanobacteria"    "Cyanobacteriia"  
#7aaa3ecae849160810cc9fedc59df606 "Bacteria" "Cyanobacteria"    "Cyanobacteriia " 
#ef9faac34d710ef9840e9a1afa1ee09f "Bacteria" "Cyanobacteria"    "Nostocophycideae"
#6f0967215a05e1eeb02ad466bea710bf "Bacteria" "Cyanobacteria"    "Nostocophycideae"
#$141f55e90b2532c91dabda3907f973e "Bacteria" "Patescibacteria " "ABY1"            
#                                 order                         family                       
#d0a102d533aa6170c2505bbef5ec7101 NA                            NA                           
#7aaa3ecae849160810cc9fedc59df606 "Cyanobacteriales "           "Chroococcidiopsaceae "      
#ef9faac34d710ef9840e9a1afa1ee09f "Nostocales"                  "Nostocaceae"                
#6f0967215a05e1eeb02ad466bea710bf "Nostocales"                  "Nostocaceae"                
#b141f55e90b2532c91dabda3907f973e "Candidatus Komeilibacteria " "Candidatus Komeilibacteria "

phy_tree(BI)
#Phylogenetic tree with 1267 tips and 1250 internal nodes.

#Tip labels:
#  d0a102d533aa6170c2505bbef5ec7101, 7aaa3ecae849160810cc9fedc59df606, ef9faac34d710ef9840e9a1afa1ee09f, 6f0967215a05e1eeb02ad466bea710bf, b141f55e90b2532c91dabda3907f973e, a8c7f7bffc30d2eebc9455dd5850ec2b, ...
#Node labels:
#  0.618, 0.335, 0.933, 0.904, 0.983, 0.853, ...

#Rooted; includes branch lengths.

taxa_names(BI) [1:10]
#[1] "d0a102d533aa6170c2505bbef5ec7101" "7aaa3ecae849160810cc9fedc59df606"
#[3] "ef9faac34d710ef9840e9a1afa1ee09f" "6f0967215a05e1eeb02ad466bea710bf"
#[5] "b141f55e90b2532c91dabda3907f973e" "a8c7f7bffc30d2eebc9455dd5850ec2b"
#[7] "3a4f17862d210ec1190e0c2cf1e6e508" "d086a6fc616a5d08fe4b650b5dab10ea"
#[9] "00f6e28d4a1265857845bcb6e81b5a1a" "2a7d4395aafcb2cdba107762c7eaae89"

###########################################################################################################################
#Plot of richness - Diagnosis

theme_set(theme_bw())
pal = "Set1"
scale_color_discrete <- function(palname = pal, ...){scale_colour_brewer(palette = palname, ...)}
scale_fill_discrete <- function(palname = pal, ...){scale_fill_brewer(palette = palname, ...)}

BI
#1267 taxa, bacteria only, do not trim more than that for alpha diversity

plot_richness(BI)   
#hard to read as it plots for all 79 samples. Can just include the alpha-diversity measures that we want
#and group samples based on Diagnosis and ProcessingMethod

plot_richness(BI, x = "Diagnosis", color = "Diagnosis", measures = c("Chao1", "Shannon"))
y = plot_richness(BI, x = "Processing_Method", color = "Diagnosis", measures = c("Chao1", "Shannon"))
y + labs(x = "Processing Method")
###########################################################################################################################
#Plot of richness - Diagnosis
#merge data based on true sample vs control sample
sampleData(BI)$Human_sample <- getVariable(BI, "Diagnosis") %in% c("ALCL", "CC", "Contralateral non-ALCL", "Cosmetic exchange", "Ruptured")

plot_richness(BI, x = "Human_sample", color = "Diagnosis", measures = c("Chao1", "Shannon"))

#merge the samples with the different diagnosis
GPst = merge_samples(BI, "Diagnosis")
#repair variables that were damaged during merge (coerced to numeric)
sample_data(GPst)$Diagnosis <- factor(sample_names(GPst))
sample_data(GPst)$Human_sample <- as.logical(sample_data(GPst)$Human_sample)

p = plot_richness(GPst, x = "Human_sample", color = "Diagnosis", measures = c("Chao1", "Shannon"))
p + geom_point(size = 4, alpha = 0.7) + labs(x = "Patient Sample")

#to remove the original smaller points on the plot
#first check which lists are present in p
p$layers
#[[1]]
#geom_point: na.rm = TRUE
#stat_identity: na.rm = TRUE
#position_identity 

#[[2]]
#mapping: ymax = ~value + se, ymin = ~value - se 
#geom_errorbar: na.rm = FALSE, orientation = NA, width = 0.1, width = 0.1, flipped_aes = FALSE
#stat_identity: na.rm = FALSE
#position_identity 

#we can see that the first layer [[1]] is the one specifying the original points, which are small. We can use negative
#indexing to pop it out, then add a new geom_point layer with larger point size
p$layers <- p$layers[-1]
p + geom_point(size = 4, alpha = 0.7) + labs(x = "Patient Sample") 

###########################################################################################################################
#Plot of richness - Processing_Method

a = plot_richness(BI, x = "Human_sample", color = "Processing_Method", measures = c("Chao1", "Shannon"))
a + labs(x = "Patient Sample")

#####################
library("vegan")
## Diagnosis
#boxplot of the number of OTUs and the Shannon entropy grouping the different sample types
plot_richness(BI, x="Diagnosis", measures=c("Observed", "Shannon")) + geom_boxplot()


# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = distance(BI, method="unifrac", weighted=F)
ordination = ordinate(BI, method="PCoA", distance=wunifrac_dist)
plot_ordination(BI, ordination, color="Diagnosis") + theme(aspect.ratio=1)

adonis(wunifrac_dist~sample_data(BI)$Diagnosis)

#Call:
#adonis(formula = wunifrac_dist ~ sample_data(BI)$Diagnosis) 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#                          Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
#sample_data(BI)$Diagnosis  6    1.8198 0.30330  1.3467 0.1009   0.01 **
#Residuals                 72   16.2152 0.22521         0.8991          
#Total                     78   18.0350                 1.0000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Processing Method
plot_richness(BI, x="Processing_Method", measures=c("Observed", "Shannon")) + geom_boxplot()


# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = distance(BI, method="unifrac", weighted=F)
ordination = ordinate(BI, method="PCoA", distance=wunifrac_dist)
plot_ordination(BI, ordination, color="Processing_Method") + theme(aspect.ratio=1)

adonis(wunifrac_dist~sample_data(BI)$Processing_Method)

#                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#sample_data(BI)$Processing_Method  1    0.5062 0.50624  2.2238 0.02807  0.007 **
#Residuals                         77   17.5288 0.22765         0.97193
#Total                             78   18.0350                 1.00000 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
