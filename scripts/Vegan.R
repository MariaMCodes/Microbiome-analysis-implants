## Convert physeq object into Vegan

BI
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1267 taxa and 79 samples ]
#sample_data() Sample Data:       [ 79 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 1267 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 1267 tips and 1250 internal nodes ]

library("vegan")


# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(BI) {
  sd <- sample_data(BI)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(BI) {
  OTU <- otu_table(BI)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}