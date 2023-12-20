setwd("../")

th <- 1000

read16s <- as.data.frame(readRDS('Table/07_prokaryote_seqtab.rds'))
readits <- as.data.frame(readRDS('Table/07_fungi_seqtab.rds')) 

th16s <- subset(read16s, rowSums(read16s) > th -1)
thits <- subset(readits, rowSums(readits) > th -1)

filt16s <- t(subset(t(th16s), rowSums(t(th16s)) > 0))
filtits <- t(subset(t(thits), rowSums(t(thits)) > 0))

taxa16S <- as.data.frame(readRDS('Table/09_prok_dominatList.rds'))
taxaITS <- as.data.frame(readRDS('Table/09_fungi_dominatList.rds'))

pro <- taxa16S[colnames(filt16s),]
fun <- taxaITS[colnames(filtits),]

write.table(pro, file='Taxa_Pro_filt.txt', quote=F, sep='\t', row.names=F)
write.table(fun, file='Taxa_Fun_filt.txt', quote=F, sep='\t', row.names=F)
