############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("../")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('tidyr', 'dplyr'))
library("igraph")
library(fungaltraits)

# -- Create directory to save
dir <- make.dir('04_Functions/05_output')

# -- Load data table
taxaITS <- readRDS('Table/09_fungi_dominatList.rds')

fngFunc <- read.csv('Table/FungalTraits 1.2_ver_16Dec_2020.csv', header=TRUE)
colnames(fngFunc)[6] <- "Genus"


############################################################################

## Fungal Trait

fnglist <- merge(taxaITS, fngFunc, by='Genus', all=FALSE, all.x=TRUE)
write.table(fnglist, file=sprintf('%s/0405_FungalTrait_all.txt', dir$tabledir), sep='\t', quote=F, row.names=F)

saveRDS(fnglist, 'Table/0405_FungalTrait_all.rds')

############################################################################