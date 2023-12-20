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

nCategory <- 6

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph'))
library(foreach)
library(ape)
library(Biostrings)

source('08_Core_microbes/functions.R')

# -- Create directory to save
dir <- make.dir('08_Core_microbes/08_output')

# -- Load data table
taxa16S <- as.data.frame(readRDS('Table/09_prok_dominatList.rds'))
taxaITS <- as.data.frame(readRDS('Table/09_fungi_dominatList.rds'))
taxa <- rbind(taxa16S, taxaITS)

dnapro <- readDNAStringSet("Table/04_OTUseq_0.98.fasta")
dnafng <- readDNAStringSet("Table/04_OTUseq_0.97.fasta")

dnapro <- data.frame(Code=names(dnapro), Sequence=as.character(dnapro))
dnafng <- data.frame(Code=names(dnafng), Sequence=as.character(dnafng))

p <- read.table('Table/07_prokaryote_annotation.txt', header=T)
f <- read.table('Table/07_fungi_annotation.txt', header=T)

rownames(p) <- p$Code
rownames(f) <- f$Code

pro <- data.frame(p, dnapro[rownames(p),])
fng <- data.frame(f, dnafng[rownames(f),])

otu.code <- data.frame(rbind(pro, fng))
rownames(otu.code) <- otu.code$ID


prefer <- foreach(i=1:nCategory) %do% {
  readRDS(sprintf('Table/08_07_preference_to_DL_Category_%s.rds', i))
}

L1 <- foreach(i=1:nCategory) %do% {
  subset(prefer[[i]], prefer[[i]]$level=='Level 1')
}

scores <- foreach(i=1:nCategory) %do% {
  tmp <- L1[[i]][,c(2, 7:9)]
  colnames(tmp) <- c('ID', 'Zvalue', 'Pvalue', 'FDR')
  rownames(tmp) <- tmp$ID
  tmp
} 

table <- foreach(i=1:nCategory) %do%{
  df <- data.frame(scores[[i]], otu.code[rownames(scores[[i]]),])
  df[,-c(5:6)]
}


foreach(i=1:nCategory) %do% {
  saveRDS(table[[i]], file=sprintf('%s/Specificity_Category_%s.rds', dir$rdsdir, i))
}

foreach(i=1:nCategory) %do% {
  write.table(table[[i]], file=sprintf('%s/Specificity_Category_%s.txt', dir$tabledir, i), quote=F, sep='\t', row.names = F)
}

######################################################
# table including all categories

name <- foreach(i=1:nCategory, .combine='c') %do% {
  rownames(table[[i]])
}

all <- unique(name)

alldf <- foreach(i=1:nCategory, .combine='cbind') %do% {
  tmp <- scores[[i]][all,][,-1]
  colnames(tmp) <- c(sprintf('Zvalue.%s', i), sprintf('Pvalue.%s', i), sprintf('FDR.%s', i))
  tmp
}

combined <- data.frame(otu.code[rownames(alldf),-10], alldf)

write.table(combined, file=sprintf('%s/Specificity_Category.txt', dir$tabledir), quote=F, sep='\t', row.names = F)
saveRDS(combined, file=sprintf('%s/Specificity_Category.rds', dir$rdsdir))

######################################################

#allcat <- na.omit(combined)



