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
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/01_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 
clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list(clr16s, clrits)
  
siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

dim(sml)
dim(dflist[[1]])
dim(dflist[[2]])

dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})

taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################

## -- Merge level
taxalevel <- 'Family'
n.core <- parallel::detectCores()

############################################################################

mergelist <- lapply(dflist, Taxa.mat, taxa=taxa, taxaLabel=taxalevel)

############################################################################
############################################################################
# euclid distance

naomit <- lapply(mergelist, function(x){ na.omit(x) })

rpca <- lapply(naomit, prcomp, scale=T)

saveRDS(rpca, 'Table/prcomp_Family.rds')

dim(naomit[[1]])
dim(naomit[[2]])

saveRDS(naomit[[1]], 'Table/comm_site10more16S_family.rds')
saveRDS(naomit[[2]], 'Table/comm_site10moreITS_family.rds')
