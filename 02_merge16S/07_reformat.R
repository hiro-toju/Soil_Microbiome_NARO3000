############################################################################
####
#### R script for Fujita (2019)
####
#### Check standard DNA quolity
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/02_merge16S")
#### 
############################################################################

## ||||||||||||||||||||||||||||||||||||||| ##
## -- Load data table
taxaprint <- readRDS("Table/06_annotation_taxonomy/taxonomy_list.rds")
seqtab <- readRDS("Table/04_seqOTUtab_0.98.rds")

## ||||||||||||||||||||||||||||||||||||||| ##
## -- Check structure
head(taxaprint)
dim(taxaprint)

dim(seqtab)
all(colnames(seqtab)==rownames(taxaprint))

## ======================================= ##
## -- Remove unidentified phylum species and remain standard DNA.
taxaBac <- taxaprint[taxaprint[,'Kingdom']%in%c('Bacteria','Archaea'),]

unident <- which(taxaBac[,'Class']=='Unidentified')
std <- grep('STD', taxaBac[,'Phylum'])

# -- Excluding Chloroplast and Mitochondria
mito = which(taxaBac == "Mitochondria", arr.ind=T)
chlo = which(taxaBac == "Chloroplast", arr.ind=T)

taxaprintIdent <- taxaBac[-c(unident, std, mito, chlo),]
dim(taxaprintIdent)

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
Bac <- taxaprintIdent[taxaprintIdent[,'Kingdom']=='Bacteria', ]
Bac[,'ID'] <- sprintf('%s_%s', unique(substr(Bac[,'Kingdom'], 1, 3)),
                  formatC(1:nrow(Bac), width = nchar(nrow(Bac)), flag = "0"))

Arc <- taxaprintIdent[taxaprintIdent[,'Kingdom']=='Archaea', ]
Arc[,'ID'] <- sprintf('%s_%s', unique(substr(Arc[,'Kingdom'], 1, 3)),
                      formatC(1:nrow(Arc), width = nchar(nrow(Arc)), flag = "0"))

BacArc <- rbind(Bac, Arc)

## ======================================= ##
## -- Remaking sequence read table
reordered <- seqtab[,rownames(BacArc)]
all(colnames(reordered)==rownames(BacArc))
colnames(reordered) <- BacArc[,'ID']

## ======================================= ##

write.table(cbind(Code=rownames(BacArc), BacArc), file='Table/07_prokaryote_annotation.txt', quote=F, sep='\t', row.names = FALSE)
saveRDS(BacArc, 'Table/07_prokaryote_annotation.rds')
saveRDS(reordered, 'Table/07_prokaryote_seqtab.rds')
saveRDS(BacArc, 'Table/prokaryote_annotation.rds')

write.csv(BacArc, 'Table/07_prokaryote_annotation.csv')
write.csv(reordered, 'Table/07_prokaryote_seqtab.csv')

## ======================================= ##
