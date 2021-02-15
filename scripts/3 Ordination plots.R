### SEQUENCING ANALYSIS USING PHYLOSEQ
## Section 3: Ordination Plots

#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")

bis
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 21 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 21 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 21 tips and 20 internal nodes ]
bis1 <- bis
bis1
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 21 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 21 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 21 tips and 20 internal nodes ]

theme_set(theme_bw())

#1) Just OTUs
BI.ord <- ordinate(bis1, "NMDS", "bray")
p1 = plot_ordination(bis1, BI.ord, type = "taxa", color = "phylum", title = "Taxa")
print(p1)

#facetting as plot p1 is hard to understand
p1 + facet_wrap(~phylum, 3)

#2) Just samples
p2 = plot_ordination(bis1, BI.ord, type = "samples", color = "Diagnosis", shape = "Human_sample")
p2 + geom_polygon(aes(fill = Diagnosis)) + geom_point(size = 4) + ggtitle("Samples")
#plot tells us nothing, dont use

#3) Biplot graphic
p3 = plot_ordination(bis1, BI.ord, type = "biplot", color = "Diagnosis", shape = "phylum", title = "biplot")
BI1.shape.names = get_taxa_unique(bis1, "phylum")
BI1.shape <- 15:(15 + length(BI1.shape.names) -1)
names(BI1.shape) <- BI1.shape.names
BI1.shape["samples"] <- 16
p3 + scale_shape_manual(values = BI1.shape)

#4) Split graphic
p4 = plot_ordination(bis1, BI.ord, type = "split", color = "phylum", shape = "Human_sample", 
                     label = "Diagnosis", title = "split")
p4

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n+1)
  hcl(h = hues, l = 65, c = 100)[1:n]
  }
color.names <- levels(p4$data$phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
p4 + scale_color_manual(values = p4cols)

###############################################################################################
##MDS ("PCoA") on Unifrac Distances
ordu = ordinate(bis1, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(bis1, ordu, color = "Diagnosis")

#to make the graph look nicer
pp = plot_ordination(bis1, ordu, color = "Diagnosis")
pp = pp + geom_point(size = 4, alpha = 0.75)
pp = pp + scale_colour_brewer(type = "qual", palette = "Set1")
pp + ggtitle("MDS/PCoA on weighted-UniFrac distance, Breast Implants")

#to remove the original smaller points on the plot as done in section 1
#first check which lists are present in p
pp$layers
#[[1]]
#geom_point: na.rm = TRUE
#stat_identity: na.rm = TRUE
#position_identity 

#[[2]]
#geom_point: na.rm = FALSE, size = 4, alpha = 0.75
#stat_identity: na.rm = FALSE
#position_identity 

#we can see that the first layer [[1]] is the one specifying the original points, which are small. We can use negative
#indexing to pop it out, then add a new geom_point layer with larger point size
pp$layers <- pp$layers[-1]
pp + geom_point(size = 4, alpha = 0.75)





