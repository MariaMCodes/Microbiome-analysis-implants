install.packages('tidyverse')
library("tidyverse")

## Counting the number of samples in each Diagnosis using plyr
install.packages('plyr')
library("plyr")

data = read.csv('data/alcl_cc2.csv')
data
count(data, 'Diagnosis')

##              Diagnosis freq
#1                   ALCL   22
#2                     CC   21
#3 Contralateral non-ALCL   15
#4      Cosmetic exchange   14
#5       Negative control    1
#6       Positive control    1
#7               Ruptured    5


## Counting the number of samples in each Diagnosis using dplyr
install.packages('dplyr')
library("dplyr")

data2 = read.csv('data/alcl_cc2.csv')

count(data2, 'Diagnosis')
count(data2, 'Processing_Method')
##  Processing_Method freq
#1                new   43
#2                old   36



