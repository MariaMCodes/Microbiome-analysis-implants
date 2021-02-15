### SEQUENCING ANALYSIS USING PHYLOSEQ
## Section 4: Introduction to Decontam

#Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")  

BiocManager::install("decontam")

#Updated the metadata file, to include a new column "Sample_or_control" 

#############################################################################################

#Check package versions
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("dplyr")
library("devtools")
library("decontam")

#Define a default theme for ggplot graphics
theme_set(theme_bw())

#import data. Note: items within data eg sample names need to match. Also had to convert data from
#biom files in csv files in conda as I had issues reading the biom files into R
otu_table = read.csv("data/table-silva_biom.csv", sep=",", row.names = 1)
otu_table

tax_table = read.csv("data/taxonomy.csv", sep=",", row.names=1, na = "")   #so R will automatically convert all the “NULL” entries in the dataset into NA
tax_table

metadata = read.csv("data/alcl_cc3.csv", sep=",", row.names=1)   #with sample or control column added to metadata
metadata
#check what structure data is in
class(metadata)

#convert to matrix
OTU <- as.matrix(otu_table)
TAX <- as.matrix(tax_table)


OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
samples = sample_data(metadata)
phy_tree = read_tree("data/rooted-tree.nwk")

#Make phyloseq object with matrix
phy <- phyloseq(OTU, TAX, samples, phy_tree)
phy
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1341 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1341 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1341 tips and 1323 internal nodes ]

#double check row names are the same - needed to add 'X' in front of sample ID and change the hiphen to a dot as biom files had this
rownames(metadata)
sample_names(phy)
#they are now the same


# 1) Filter out everything but the kingdom of Bacteria, since it is not our intent to amplify anything else
phy_sub <- subset_taxa(phy, kingdom == "Bacteria")
phy_sub
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]
#lost 74 taxa from the OTU table, meaning out amplification seems to have worked well to only target bacteria

#rename data to track easily (DD = decontam data)
DD <- phy_sub
DD 

# 2) Inspect library sizes
#first look at the library sizes (ie the number of reads) in each sample, as a function of whether that sample was a true
#positive sample or a negative control

#see depths
df <- as.data.frame(sample_data(DD))     #put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(DD)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) + geom_point()


###############################################################################################################
##Removing contaminants
# 3) Identify contaminants - Prevalence
sample_data(DD)$is.neg <- sample_data(DD)$Sample_or_Control == "Control sample"
contamdf.prev <- isContaminant(DD, method = "prevalence", neg = "is.neg")
table(contamdf.prev$contaminant)

#FALSE  TRUE 
#1252    15  

which(contamdf.prev$contaminant)
#[1]  225  391  432  471  512  576  638  699  711  733  812  974  995 1013 1224

#prevalence-based contaminant identification has identified 15 contaminants in this dataset. This is the default method
#of classification threshold

#try using the more aggressive classification threshold
contamdf.prev05 <- isContaminant(DD, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev05$contaminant)
#FALSE  TRUE 
#1244    23 
which(contamdf.prev05$contaminant)
#[1]  150  225  390  391  404  420  432  471  502  512  576  638  699  711  733  812  974
#[18]  995 1013 1114 1179 1223 1224

#try using less aggressive classification thresholds
#0.4
contamdf.prev04 <- isContaminant(DD, method = "prevalence", neg = "is.neg", threshold = 0.4)
table(contamdf.prev04$contaminant)
#FALSE  TRUE 
#1246    21  

#0.3
contamdf.prev03 <- isContaminant(DD, method = "prevalence", neg = "is.neg", threshold = 0.3)
table(contamdf.prev03$contaminant)
#FALSE  TRUE 
#1248    19 

#0.2
contamdf.prev02 <- isContaminant(DD, method = "prevalence", neg = "is.neg", threshold = 0.2)
table(contamdf.prev02$contaminant)
#FALSE  TRUE 
#1250    17  

#0.25
contamdf.prev025 <- isContaminant(DD, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev025$contaminant)

#################################################################################################

#look at the number of times several of these taxa were observed in negative controls and positive samples
#make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(DD, function(abund) 1*(abund > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True sample", ps.pa)

#Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg), 
                    contaminant = contamdf.prev03$contaminant)                      #this is a graph of the more aggressive classification threshold
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Control)") + ylab("Prevalence (True Samples)") 

# 4) Now that we have identified likely contaminants, let's remove them from the phyloseq object

DD
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]

#using the default classification threshold
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, DD)   
#need to use original DD data otherwise the ps.pa data is not in the OTU_table format since you transformed it ie function(abund) 1*(abund > 0)
ps.noncontam
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1252 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1252 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1252 tips and 1237 internal nodes ]

#15 taxa removed - this is contaminant-filtered data

#using the more aggressive classification threshold
ps.noncontam2 <- prune_taxa(!contamdf.prev05$contaminant, DD)  
ps.noncontam2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1244 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1244 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1244 tips and 1232 internal nodes ]

#23 taxa removed - this is contaminant-filtered data
#herein using ps.noncontam2

#0.4
ps.noncontam_4 <- prune_taxa(!contamdf.prev04$contaminant, DD)  
ps.noncontam_4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1246 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1246 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1246 tips and 1233 internal nodes ]

#using the less aggressive classification threshold (0.3)
ps.noncontam3 <- prune_taxa(!contamdf.prev03$contaminant, DD)  
ps.noncontam3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1248 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1248 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1248 tips and 1234 internal nodes ]

#0.2
ps.noncontam4 <- prune_taxa(!contamdf.prev02$contaminant, DD)  
ps.noncontam4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1250 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1250 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1250 tips and 1235 internal nodes ]

#0.25
ps.noncontam25 <- prune_taxa(!contamdf.prev025$contaminant, DD)  
ps.noncontam25
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1244 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1244 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1244 tips and 1232 internal nodes ]

###########################################################################################
##Preprocessing steps
#a) filtering
DDr = transform_sample_counts(ps.noncontam2, function(x) x/sum(x))
DDfr = filter_taxa(DDr, function(x) mean(x) > 1e-5, TRUE)
DDfr
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1162 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1162 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1162 tips and 1150 internal nodes ]

DDr2 = transform_sample_counts(ps.noncontam2, function(x) x/sum(x))
DDfr2 = filter_taxa(DDr2, function(x) mean(x) > 1e-3, TRUE)
DDfr2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 158 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 158 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 158 tips and 156 internal nodes ]

#b) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
DD2 = filter_taxa(ps.noncontam2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
DD2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 13 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 13 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 13 tips and 12 internal nodes ]

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(DD2))
standf = function(x, t = total) round(t*(x/sum(x)))
dds = transform_sample_counts(DD2, standf)
dds

#d) filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
DDfr3 = filter_taxa(dds, function(x) sd(x)/mean(x) > 3.0, TRUE)
DDfr3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1 taxa by 6 taxonomic ranks ]
##will remove this step since only one taxa left 

##from step c)
#plot the only 13 taxa i.e. everything
plot_bar(dds, "Diagnosis", fill = "kingdom")
plot_bar(dds, "Diagnosis", fill = "phylum")
plot_bar(dds, "Diagnosis", fill = "class")
plot_bar(dds, "Diagnosis", fill = "order")
plot_bar(dds, "Diagnosis", fill = "family")
plot_bar(dds, "Diagnosis", fill = "genus")


#################################################################################################
#since pseudomonas disappeared from positive control we will try the default classification threshold

ps.noncontam
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1252 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1252 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1252 tips and 1237 internal nodes ]


##Preprocessing steps
#b) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
DD3 = filter_taxa(ps.noncontam, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
DD3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 21 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 21 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 21 tips and 20 internal nodes ]

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(DD3))
standf = function(x, t = total) round(t*(x/sum(x)))
dds2 = transform_sample_counts(DD3, standf)
dds2

#plot the only 13 taxa i.e. everything
plot_bar(dds2, "Diagnosis", fill = "kingdom")
plot_bar(dds2, "Diagnosis", fill = "phylum")
plot_bar(dds2, "Diagnosis", fill = "class")
plot_bar(dds2, "Diagnosis", fill = "order")
plot_bar(dds2, "Diagnosis", fill = "family")
plot_bar(dds2, "Diagnosis", fill = "genus")

#################################################################################################
#try the less aggressive classification threshold 0.3

ps.noncontam3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1248 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1248 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1248 tips and 1233 internal nodes ]


##Preprocessing steps
#b) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
DD4 = filter_taxa(ps.noncontam3, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
DD4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 17 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 17 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 17 tips and 16 internal nodes ]

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(DD4))
standf = function(x, t = total) round(t*(x/sum(x)))
dds3 = transform_sample_counts(DD4, standf)
dds3

#plot the only 13 taxa i.e. everything
plot_bar(dds3, "Diagnosis", fill = "kingdom")
plot_bar(dds3, "Diagnosis", fill = "phylum")
plot_bar(dds3, "Diagnosis", fill = "class")
plot_bar(dds3, "Diagnosis", fill = "order")
plot_bar(dds3, "Diagnosis", fill = "family")
plot_bar(dds3, "Diagnosis", fill = "genus")

#################################################################################################
#try the less aggressive classification threshold 0.2

ps.noncontam4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1250 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 1250 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1250 tips and 1235 internal nodes ]


##Preprocessing steps
#b) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
DD5 = filter_taxa(ps.noncontam4, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
DD5
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 19 taxa by 6 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 19 tips and 18 internal nodes ]

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(DD5))
standf = function(x, t = total) round(t*(x/sum(x)))
dds4 = transform_sample_counts(DD5, standf)
dds4

#plot the only 13 taxa i.e. everything
plot_bar(dds4, "Diagnosis", fill = "kingdom")
plot_bar(dds4, "Diagnosis", fill = "phylum")
plot_bar(dds4, "Diagnosis", fill = "class")
plot_bar(dds4, "Diagnosis", fill = "order")
plot_bar(dds4, "Diagnosis", fill = "family")
plot_bar(dds4, "Diagnosis", fill = "genus")

####################################################################################################
#try the less aggressive classification threshold 0.4

ps.noncontam_4

##Preprocessing steps
#b) remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean and trivially large CV
DD_4 = filter_taxa(ps.noncontam_4, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
DD_4
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 15 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 4 sample variables ]
#tax_table()   Taxonomy Table:    [ 15 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 15 tips and 14 internal nodes ]]

#c) standardise abundances to the median sequencing depth
total = median(sample_sums(DD_4))
standf = function(x, t = total) round(t*(x/sum(x)))
dds_4 = transform_sample_counts(DD_4, standf)
dds_4


