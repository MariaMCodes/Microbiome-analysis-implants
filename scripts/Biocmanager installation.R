### BiocManager

## Introduction
#Use the BiocManager package to install and manage packages from the Bioconductor project
#for the statistical analysis and comprehension of high-throughput genomic data.

#Current Bioconductor packages are available on a 'release' version intended for every-day use, 
#and a 'devel' version where new features are introduced. A new release version is created every 
#six months. Using the BiocManager package helps users install packages from the same release.


## Installing BiocManager
#Use standard R installation procedures to install the BiocManager package. This command is requried 
#only once per R installation.
chooseCRANmirror()
install.packages("BiocManager")
2

## Installing Bioconductor, CRAN, or github packages
#Install Bioconductor (or CRAN) packages with
BiocManager::install(c("GenomicRanges", "Organism.dplyr"))

#Installed packages can be updated to their current version with
BiocManager::install()


## Version and validity of installations
#Use version() to discover the version of Bioconductor currently in use.
BiocManager::version()

#Bioconductor packages work best when they are all from the same release. 
#Use valid() to identify packages that are out-of-date or from unexpected versions.
BiocManager::valid()

#valid() returns an object that can be queried for detailed information about invalid packages


## Available packages
#Packages available for your version of Bioconductor can be discovered with available(); the first 
#argument can be used to filter package names based on a regular expression, e.g., 'BSgenome' package 
#available for Homo sapiens
avail <- BiocManager::available()
length(avail)                               # all CRAN & Bioconductor packages
BiocManager::available("BSgenome.Hsapiens") # BSgenome.Hsapiens.* packages


