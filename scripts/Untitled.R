### Importing phyloseq Data

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

#Define a default theme for ggplot graphics.
theme_set(theme_bw())

## phyloseq-ize Data already in R
#Any data already in an R session can be annoated/coerced to be recognized by phyloseq’s functions and methods. 
#This is important, because there are lots of ways you might receive data related to a microbiome project, and 
#not all of these will come from a popular server or workflow that is already supported by a phyloseq import function. 
#There are, however, lots of general-purpose tools available for reading any common file format into R. We also want to 
#encourage users to create and share their own import code for special/new data formats as they arive.

#For these reasons especially, phyloseq provides tools for constructing phyloseq component data, and the experiment-level 
#multi-component data object, the phyloseq-class. These are the same functions used internally by the currently available importers.
?phyloseq

## Build phyloseq-class objects from their components.
#phyloseq() is a constructor method, This is the main method suggested for constructing an experiment-level (phyloseq-class) object from 
#its component data (component data classes: otu_table-class, sample_data-class, taxonomyTable-class, phylo-class).

#Usage
#phyloseq(...)
#Arguments
#...	- One or more component objects among the set of classes defined by the phyloseq package, as well as phylo-class (defined by the ape-package). 
#Each argument should be a different class. For combining multiple components of the same class, or multiple phyloseq-class objects, use the merge_phyloseq 
#function. Unlike in earlier versions, the arguments to phyloseq do not need to be named, and the order of the arguments does not matter.

#Value
#The class of the returned object depends on the argument class(es). For an experiment-level object, two or more component data objects must be provided. 
#Otherwise, if a single component-class is provided, it is simply returned as-is. The order of arguments does not matter.

#Examples
#data(esophagus)
#x1 = phyloseq(otu_table(esophagus), phy_tree(esophagus))
identical(x1, esophagus)
# # data(GlobalPatterns)
# # GP <- GlobalPatterns
# # phyloseq(sample_data(GP), otu_table(GP))
# # phyloseq(otu_table(GP), phy_tree(GP))
# # phyloseq(tax_table(GP), otu_table(GP))
# # phyloseq(phy_tree(GP), otu_table(GP), sample_data(GP))
# # phyloseq(otu_table(GP), tax_table(GP), sample_data(GP))
# # phyloseq(otu_table(GP), phy_tree(GP), tax_table(GP), sample_data(GP))

?otu_table
## Build or access the otu_table.
#This is the suggested method for both constructing and accessing Operational Taxonomic Unit (OTU) abundance (otu_table-class) objects. When the first 
#argument is a matrix, otu_table() will attempt to create and return an otu_table-class object, which further depends on whether or not taxa_are_rows 
#is provided as an additional argument. Alternatively, if the first argument is an experiment-level (phyloseq-class) object, then the corresponding 
#otu_table is returned.

#Usage
#otu_table(object, taxa_are_rows, errorIfNULL=TRUE)

## S4 method for signature 'phyloseq'
#otu_table(object, errorIfNULL = TRUE)

## S4 method for signature 'otu_table'
#otu_table(object, errorIfNULL = TRUE)

## S4 method for signature 'matrix'
#otu_table(object, taxa_are_rows)

## S4 method for signature 'data.frame'
#otu_table(object, taxa_are_rows)

## S4 method for signature 'ANY'
#otu_table(object, errorIfNULL = TRUE)

#Arguments
#object	- (Required). An integer matrix, otu_table-class, or phyloseq-class.
#taxa_are_rows - (Conditionally optional). Logical; of length 1. Ignored unless object is a matrix, in which case it is is required.
#errorIfNULL - (Optional). Logical. Should the accessor stop with an error if the slot is empty (NULL)? Default TRUE. Ignored if object 
#argument is a matrix (constructor invoked instead).

#Value
#An otu_table-class object.

#See Also
#phy_tree, sample_data, tax_table phyloseq, merge_phyloseq

#Examples
# data(GlobalPatterns)
# otu_table(GlobalPatterns)

?sample_data
## Build or access sample_data.
#This is the suggested method for both constructing and accessing a table of sample-level variables (sample_data-class), which in the 
#phyloseq-package is represented as a special extension of the data.frame-class. When the argument is a data.frame, sample_data will 
#create a sample_data-class object. In this case, the rows should be named to match the sample_names of the other objects to which it 
#will ultimately be paired. Alternatively, if the first argument is an experiment-level (phyloseq-class) object, then the corresponding 
#sample_data is returned. Like other accessors (see See Also, below), the default behavior of this method is to stop with an error if object 
#is a phyloseq-class but does not contain a sample_data.

#Usage
#sample_data(object, errorIfNULL=TRUE)

## S4 method for signature 'ANY'
#sample_data(object, errorIfNULL = TRUE)

## S4 method for signature 'data.frame'
#sample_data(object)

#Arguments
#object - (Required). A data.frame-class, or a phyloseq-class object.
#errorIfNULL - (Optional). Logical. Should the accessor stop with an error if the slot is empty (NULL)? Default TRUE.

#Value
#A sample_data-class object representing the sample variates of an experiment.

#See Also
#phy_tree, tax_table, otu_table phyloseq, merge_phyloseq

Examples
#data(soilrep)
#head(sample_data(soilrep))

?tax_table
## Build or access the taxonomyTable.
#This is the suggested method for both constructing and accessing a table of taxonomic names, organized with ranks as columns (taxonomyTable-class). When the argument is a character matrix, tax_table() will create and return a taxonomyTable-class object. In this case, the rows should be named to match the species.names of the other objects to which it will ultimately be paired. Alternatively, if the first argument is an experiment-level (phyloseq-class) object, then the corresponding taxonomyTable is returned. Like other accessors (see See Also, below), the default behavior of this method is to stop with an error if object is a phyloseq-class but does not contain a taxonomyTable.

#Usage
#tax_table(object, errorIfNULL=TRUE)

## S4 method for signature 'ANY'
#tax_table(object, errorIfNULL = TRUE)

## S4 method for signature 'matrix'
#tax_table(object)

## S4 method for signature 'data.frame'
#tax_table(object)

#Arguments
#object	- An object among the set of classes defined by the phyloseq package that contain taxonomyTable.
#errorIfNULL - (Optional). Logical. Should the accessor stop with an error if the slot is empty (NULL)? Default TRUE.

#Value
#A taxonomyTable-class object. It is either grabbed from the relevant slot if object is complex, or built anew if object 
#is a character matrix representing the taxonomic classification of species in the experiment.

#See Also
#phy_tree, sample_data, otu_table phyloseq, merge_phyloseq

#Examples
# tax1 <- tax_table(matrix("abc", 30, 8))
# data(GlobalPatterns)
# tax_table(GlobalPatterns)

