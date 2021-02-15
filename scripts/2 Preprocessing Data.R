### SEQUENCING ANALYSIS USING PHYLOSEQ
## SECTION 2: Pre-Processing Data

#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

BI

#####################################################################################################
#Test - filter low-occurence, poorly represented OTUs from the data, because they are essentially noise
#variables. Remove OTUs that do not shoe appear moer than 5 times in more than half the samples
wh0 = genefilter_sample(BI, filterfun_sample(function(x) x > 5), A = 0.5*nsamples(BI))
BI1 = prune_taxa(wh0, BI)
BI1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 4 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4 tips and 3 internal nodes ]
#see what they are
tax_table(BI) [1:10, 1:6] #got some weird genus eg Parcubacteria so not a good filter
#####################################################################################################

#1) Filtering data
#the data is first transformed to relative abundance, creating the new BIr object, which is then filtered
#such that only OTUs with a mean greater than 10^-5 are kept
BIr = transform_sample_counts(BI, function(x) x / sum(x))  # transformation
BIr = filter_taxa(BIr, function(x) mean(x) > 1e-5, TRUE)   # filtering
BIr
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1122 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1122 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1122 tips and 1105 internal nodes ]

#this results in a subsetted object, BIr, containing just 1122 taxa of the original 1267 OTUs

#2) Subsetting data 
#we can assign to BI.chl the subset of the BI dataset that are part of the Chlamydiae phylum, and then remove 
#samples with less than 20 total reads
BI.chl = subset_taxa(BI, phylum == "Chlamydiae")
BI.chl = prune_samples(sample_sums(BI.chl) >= 20, BI.chl)
BI.chl
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4 taxa and 3 samples ]
#sample_data() Sample Data:       [ 3 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 4 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4 tips and 3 internal nodes ]

#3) Merging data
#to merge all OTUs in the Chlamydiae-only dataset (only 4 OTUs). Can also add [1:4]
BI.chl.merged = merge_taxa(BI.chl, taxa_names(BI.chl))   

#4) Agglomeration functions
#to agglomerate the "Chlamydiae-only" dataset (called bichla) at the taxonomic rank of Family, and create
#an annotated tree as a result
bichla = tax_glom(BI.chl, "family")
plot_tree(bichla, color = "Diagnosis", shape = "class", size = "Abundance")
plot_tree(bichla, color = "Diagnosis", shape = "order", size = "Abundance")
plot_tree(bichla, color = "Diagnosis", shape = "family", size = "Abundance")
plot_tree(bichla, color = "Diagnosis", shape = "genus", size = "Abundance")


#5) Transforming abundance values
#to transform BI.chl abundance counts to fractional abundance
transform_sample_counts(BI.chl, function(OTU) OTU/sum(OTU))

###############################################################################################################################################
#6) Remaining preprocessing steps
#a) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
BI2 = filter_taxa(BI, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
BI2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 21 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 21 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 21 tips and 20 internal nodes ]

#did this in previous section (Section 1) for alpha diversity measures
#b) define a human versus non-human categorical variable, and add this new variable to sample data
#sample_data(BI2)$human = factor(get_variable(BI2, "Diagnosis") %in% c("ALCL", "Contralateral non-ALCL",
#                                                                    "CC", "Ruptured", "Cosmetic exchange"))

sample_variables(BI2)
#[1] "Diagnosis"        "ProcessingMethod" "Human_sample"       #now "human" has been added"

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(BI2))
standf = function(x, t = total) round(t * (x/sum(x)))
bis = transform_sample_counts(BI2, standf)
bis

#d) filter the taxa using a cutoff of 3.0 for the Coefficient of Variation, which results in only 2 taxa left
bisf = filter_taxa(bis, function(x) sd(x)/mean(x) > 3.0, TRUE)
bisf

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 2 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 2 tips and 1 internal nodes ]

#check table with only 2 taxa
tax_table(bisf) [1:2]
#Taxonomy Table:     [2 taxa by 6 taxonomic ranks]:
#                                 kingdom    phylum           class                  order             
#352d6677c8563549ac3279667bb0505a "Bacteria" "Proteobacteria" "Gammaproteobacteria " "Pseudomonadales "
#7c82e9683757ca5c546929207e33e95a "Bacteria" "Firmicutes"     "Bacilli"              "Bacillales "     
#                                 family           genus          
#352d6677c8563549ac3279667bb0505a "Moraxellaceae " "Enhydrobacter"
#7c82e9683757ca5c546929207e33e95a "Bacillaceae"    "Geobacillus"  

#e) subset the data to "Proteobacteria"
bisfp = subset_taxa(bisf, phylum == "Proteobacteria")

#7) Graphic summary
title = "Proteobacteria-only"
plot_bar(bisfp, "Diagnosis", "Abundance", title = title)
plot_bar(bisfp, "Diagnosis", "Abundance", "genus", title = title)

#subset the data to "Firmicutes"
bisff = subset_taxa(bisf, phylum == "Firmicutes")
title = "Firmicutes-only"
plot_bar(bisff, "Diagnosis", "Abundance", title = title)
plot_bar(bisff, "Diagnosis", "Abundance", "family", title = title)

#plot the only 2 taxa i.e. everything
plot_bar(bisf, "Diagnosis", fill = "genus")

#####################################################################################################
####  IF DON'T DO STEP d)  ####
#working on only bis data
#check table with 21 taxa
tax_table(bis) [1:21]

#e) subset the data to "Proteobacteria" 
bisp = subset_taxa(bis, phylum == "Proteobacteria")

#7) Graphic summary
title = "Proteobacteria-only"
plot_bar(bisp, "Diagnosis", "Abundance", title = title)
plot_bar(bisp, "Diagnosis", "Abundance", "genus", title = title)

#subset the data to "Firmicutes"
bisf = subset_taxa(bis, phylum == "Firmicutes")
title = "Firmicutes-only"
plot_bar(bisf, "Diagnosis", "Abundance", title = title)
plot_bar(bisf, "Diagnosis", "Abundance", "family", title = title)

#plot the only 21 taxa i.e. everything
plot_bar(bis, "Diagnosis", fill = "kingdom")
plot_bar(bis, "Diagnosis", fill = "phylum")
plot_bar(bis, "Diagnosis", fill = "class")
plot_bar(bis, "Diagnosis", fill = "order")
plot_bar(bis, "Diagnosis", fill = "family")
plot_bar(bis, "Diagnosis", fill = "genus")
plot_bar(bis, "Abundance", fill = "genus", facet_grid = "Diagnosis~.")

##################################################################################################
##how to expand color palettes in R
library("RColorBrewer")
#Define the number of colors you want
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
#Create bar plot with 15 colors
#Use scale_fill_manual
plot_bar(bis, "Diagnosis", fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + theme_bw() +
  theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#Plot of Gram neg vs Gram pos
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
plot_bar(bis, "Diagnosis", fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + theme_bw() +
  theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#more sophisticated organization using facets
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
plot_bar(bis, "phylum", fill = "genus", facet_grid = ~Diagnosis)

plot_bar(bis, "Diagnosis", fill = "genus", facet_grid = ~Processing_Method)

##################################################################################################
##To remove the black lines from barplot
#look at the code for plot_bar
plot_bar
#function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
#title = NULL, facet_grid = NULL) 
#{
#  mdf = psmelt(physeq)
#  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
#  p = p + geom_bar(stat = "identity", position = "stack", color = "black")
#  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
#  if (!is.null(facet_grid)) {
#    p <- p + facet_grid(facet_grid)
#  }
#  if (!is.null(title)) {
#    p <- p + ggtitle(title)
#  }
#  return(p)
#}
#<bytecode: 0x7fd476596f20>
#  <environment: namespace:phyloseq>

#As you can see, the function includes this statement:
#geom_bar(stat = "identity", position = "stack", color = "black")

#The color="black" argument is what causes the black lines. This is a pretty basic 
#bar plot and you can just create your own function based on this code:

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
          title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

#Notice that the only change is that color="black" has been removed. Now run my_plot_bar 
#instead of plot_bar and get a bar plot without the black lines.

#Grouped by diagnosis
library("RColorBrewer")
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(bis, "Diagnosis", fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#all samples
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(bis, fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

#Grouped by Gram positive or negative
library("RColorBrewer")
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
my_plot_bar(bis, "Diagnosis", fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#all samples
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
my_plot_bar(bis, fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#Facetted by diagnosis
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(bis, "genus", fill = "genus", facet_grid = "Diagnosis~.") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")
#to have no color
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(bis, "genus", facet_grid = "Diagnosis~.") +
  geom_bar(stat = "identity", position = "stack")  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Genus")

#Grouped by Processing Method
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(bis, "Diagnosis", fill = "genus", facet_grid = ~Processing_Method) + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")


#useful website on tips to edit legend https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/

