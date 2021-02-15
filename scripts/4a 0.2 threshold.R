### SEQUENCING ANALYSIS USING PHYLOSEQ
## Dection 4a: Analysis on default classification method


library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")
library("decontam")


dds4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 19 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 19 tips and 18 internal nodes ]

theme_set(theme_bw())


##MDS ("PCoA") on Unifrac Distances
ordu = ordinate(dds4, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(dds4, ordu, color = "Diagnosis")

#to make the graph look nicer
pp = plot_ordination(dds4, ordu, color = "Diagnosis")
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


###################################################################################################

##how to expand color palettes in R
library("RColorBrewer")
#Define the number of colors you want
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
#Create bar plot with 15 colors
#Use scale_fill_manual
plot_bar(dds4, "Diagnosis", fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + theme_bw() +
  theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#Plot of Gram neg vs Gram pos
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
plot_bar(dds4, "Diagnosis", fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + theme_bw() +
  theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#more sophisticated organization using facets
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
plot_bar(dds4, "phylum", fill = "genus", facet_grid = ~Diagnosis)

plot_bar(dds4, "Diagnosis", fill = "genus", facet_grid = ~Processing_Method)

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
my_plot_bar(dds4, "Diagnosis", fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#all samples
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(dds4, fill = "genus") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#Grouped by Gram positive or negative
library("RColorBrewer")
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
my_plot_bar(dds4, "Diagnosis", fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#all samples
nb.rows <- 2
mycolors <- colorRampPalette(brewer.pal(4, "Set3")) (nb.rows)
my_plot_bar(dds4, fill = "Gram_Pos_or_Neg") + scale_fill_manual(values = mycolors, name = "") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")

#Facetted by diagnosis
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(dds4, "genus", fill = "genus", facet_grid = "Diagnosis~.") + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")
#to have no color
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(dds4, "genus", facet_grid = "Diagnosis~.") +
  geom_bar(stat = "identity", position = "stack")  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Genus")

#Grouped by Processing Method
nb.rows <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3")) (nb.rows)
my_plot_bar(dds4, "Diagnosis", fill = "genus", facet_grid = ~Processing_Method) + scale_fill_manual(values = mycolors, name = "Genus") + 
  theme_bw() +  theme(legend.position = "bottom", legend.title.align = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "")



