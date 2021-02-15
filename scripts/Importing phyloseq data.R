

#We’ll create the example vanilla R tables using base R code. No packages required yet.
#Create a pretend OTU table that you read from a file, called otumat
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat
#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#[1,]    44   63   35   23   97    8   65   11   37    54
#[2,]    86   64   18   17   25   40   20   22   71    65
#[3,]    67   41   54   83   68   60   98   94   75    49
#[4,]    41   73   20    9   98   73   92   70   17    93
#[5,]    42   57   46   99   91   55    2   34   25    87
#[6,]     9   10   60   58   84   64   49   97   27    23
#[7,]    98   88   83   33   29   67   49   44   33    78
#[8,]    65   75   87   57   98   91   49   86   74    78
#[9,]    65   14   90   96    4   49    3   61   34    36
#[10,]   99   18   34   97   18   37   71   87   20    33


#It needs sample names and OTU names, the index names of the your own matrix might already have this.
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat
#Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8 Sample9 Sample10
#OTU1       44      63      35      23      97       8      65      11      37       54
#OTU2       86      64      18      17      25      40      20      22      71       65
#OTU3       67      41      54      83      68      60      98      94      75       49
#OTU4       41      73      20       9      98      73      92      70      17       93
#OTU5       42      57      46      99      91      55       2      34      25       87
#OTU6        9      10      60      58      84      64      49      97      27       23
#OTU7       98      88      83      33      29      67      49      44      33       78
#OTU8       65      75      87      57      98      91      49      86      74       78
#OTU9       65      14      90      96       4      49       3      61      34       36
#OTU10      99      18      34      97      18      37      71      87      20       33

#Now we need a pretend taxonomy table
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat
#Domain Phylum Class Order Family Genus Species
#OTU1  "y"    "g"    "c"   "b"   "v"    "u"   "p"    
#OTU2  "v"    "m"    "t"   "j"   "c"    "y"   "x"    
#OTU3  "x"    "j"    "n"   "m"   "s"    "i"   "z"    
#OTU4  "g"    "e"    "r"   "t"   "w"    "a"   "j"    
#OTU5  "r"    "x"    "x"   "g"   "g"    "f"   "x"    
#OTU6  "d"    "u"    "f"   "p"   "v"    "o"   "p"    
#OTU7  "j"    "g"    "h"   "x"   "l"    "g"   "t"    
#OTU8  "y"    "n"    "b"   "k"   "o"    "s"   "g"    
#OTU9  "i"    "q"    "v"   "y"   "z"    "s"   "g"    
#OTU10 "y"    "g"    "s"   "f"   "p"    "s"   "v"  

class(otumat)

class(taxmat)


#Note how these are just vanilla R matrices. Now let’s tell phyloseq how to combine them into a 
#phyloseq object. In the previous lines, we didn’t even need to have phyloseq loaded yet. Now we do.
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
#OTU Table:          [10 taxa and 10 samples]
#                    taxa are rows
#      Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8 Sample9 Sample10
#OTU1       44      63      35      23      97       8      65      11      37       54
#OTU2       86      64      18      17      25      40      20      22      71       65
#OTU3       67      41      54      83      68      60      98      94      75       49
#OTU4       41      73      20       9      98      73      92      70      17       93
#OTU5       42      57      46      99      91      55       2      34      25       87
#OTU6        9      10      60      58      84      64      49      97      27       23
#OTU7       98      88      83      33      29      67      49      44      33       78
#OTU8       65      75      87      57      98      91      49      86      74       78
#OTU9       65      14      90      96       4      49       3      61      34       36
#OTU10      99      18      34      97      18      37      71      87      20       33

TAX
#Taxonomy Table:     [10 taxa by 7 taxonomic ranks]:
#      Domain Phylum Class Order Family Genus Species
#OTU1  "y"    "g"    "c"   "b"   "v"    "u"   "p"    
#OTU2  "v"    "m"    "t"   "j"   "c"    "y"   "x"    
#OTU3  "x"    "j"    "n"   "m"   "s"    "i"   "z"    
#OTU4  "g"    "e"    "r"   "t"   "w"    "a"   "j"    
#OTU5  "r"    "x"    "x"   "g"   "g"    "f"   "x"    
#OTU6  "d"    "u"    "f"   "p"   "v"    "o"   "p"    
#OTU7  "j"    "g"    "h"   "x"   "l"    "g"   "t"    
#OTU8  "y"    "n"    "b"   "k"   "o"    "s"   "g"    
#OTU9  "i"    "q"    "v"   "y"   "z"    "s"   "g"    
#OTU10 "y"    "g"    "s"   "f"   "p"    "s"   "v"  

physeq = phyloseq(OTU, TAX)
physeq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
#tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]

plot_bar(physeq, fill = "Family")

