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
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph', 'ggnetwork', 'graphlayouts', 'RColorBrewer'))
library(tidyverse)
library(ape)
library(Biostrings)

# -- Create directory to save
dir <- make.dir('08_Core_microbes/06_output')

# -- Load data table

x <- read.csv('Table/0804_networkSummary.csv', header=T, row.name=1)

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

combined <- left_join(x[, -c(7:10, 13:19, 23:25)], otu.code[, -10], by='ID')
combined$betweenness <- round(combined$betweenness, 5)

saveRDS(combined, 'Table/0806_NetworkOTUs_summary.rds')
write.table(combined, file=sprintf('%s/0806_NetworkOTUs_summary.txt', dir$tabledir), sep='\t', quote=F, row.names=F)

