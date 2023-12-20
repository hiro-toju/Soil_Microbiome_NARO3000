############################################################################
####
#### R script for Fujita (2019)
####
#### Check standard DNA quolity
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220320/02_mergeITS")
#### 
############################################################################

## ||||||||||||||||||||||||||||||||||||||| ##
## -- Load data table
taxaprint <- readRDS("Table/06_annotation_taxonomy/taxonomy_list.rds")
seqtab <- readRDS("Table/04_seqOTUtab_0.97.rds")

rownames(seqtab) <- gsub("S_0", "S_", rownames(seqtab))

## ||||||||||||||||||||||||||||||||||||||| ##
## -- Check structure
head(taxaprint)
dim(taxaprint)

head(seqtab[,1:10])
dim(seqtab)
all(colnames(seqtab)==rownames(taxaprint))

## ======================================= ##
## -- Compile taxaprint
taxaFng <- taxaprint[taxaprint[,'Kingdom']=='Fungi', c('Kingdom', 'Phylum','Class','Order','Family','Genus','Species')]
colnames(taxaFng) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
taxaFng[taxaFng==''] <- 'Unidentified'

## ======================================= ##
## -- Remove unidentified phylum species and remain standard DNA.
unident <- which(taxaFng[,'Class']=='Unidentified')
std <- grep('STD', taxaFng[,'Phylum'])
taxaprintIdent <- taxaFng[-c(unident, std),]

## -- Showing no effectof unident phyum to down stream 
rel <- seqtab/rowSums(seqtab)

dir.create('Figures')
pdf('Figures/notIdent_phyum_abundance.pdf')
plot(apply(rel[,c(unident, std)], 2, max, na.rm=TRUE),
     ylab='Maximum relative abundance')
plot(colSums(rel[,c(unident, std)]>0, na.rm = TRUE),
     ylab='Number of samples which unident phyum appeared in')
dev.off()


## ======================================= ##
## -- Renaming OTU IDs.
Fng <- taxaprintIdent[taxaprintIdent[,'Kingdom']=='Fungi', ]
FngID <- sprintf('%s_%s', unique(substr(Fng[,'Kingdom'], 1, 3)),
                  formatC(1:nrow(Fng), width = nchar(nrow(Fng)), flag = "0"))
Fng <- cbind(ID=FngID, Fng)

## ======================================= ##
## -- Remaking sequence read table
reordered <- seqtab[,rownames(Fng)]
all(colnames(reordered)==rownames(Fng))
colnames(reordered) <- Fng[,'ID']

## ======================================= ##

write.table(cbind(Code=rownames(Fng), Fng), file='Table/07_fungi_annotation.txt', quote=F, sep='\t', row.names = FALSE)
saveRDS(Fng, 'Table/07_fungi_annotation.rds')
saveRDS(reordered, 'Table/07_fungi_seqtab.rds')
saveRDS(reordered, 'Table/fungi_seqtab.rds')

write.csv(Fng, 'Table/07_fungi_annotation.csv')
write.csv(reordered, 'Table/07_fungi_seqtab.csv')

## ======================================= ##
