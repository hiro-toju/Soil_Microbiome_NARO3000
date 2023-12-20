##############################################
##
## 2022. 2. 14
## Extract STDsequence
##
##############################################

library(seqinr)
library(dada2)

blast <- read.table("Table/03_blastResult.txt", stringsAsFactors=FALSE)
seq <- read.fasta("Table/02_16Smerge_clus_seq.fasta")
seqtab <- readRDS("Table/02_16Smerge_seqtab.rds")

# -- Reference database
reference="Table/STD_prok.fasta"

##############################################
# -- Data compile to other taxonomy level matrix

Taxa.mat <- function(x, y, taxaLabel, func=function(x){sum(x)}){
  
  colnames(x) <- y[colnames(x), taxaLabel]
  
  summary <- do.call(cbind,
                     lapply(unique(colnames(x)), 
                            function(a){ num <- which(colnames(x)==a)
                            apply(as.matrix(x[,num]), 1, func)}) )
  
  colnames(summary) <- unique(colnames(x))
  summary
}

##############################################
alignlen <- blast[which(blast[,4] > 200 & blast[,3] > 90),]

STDseq <- seq[names(seq) %in% unique(alignlen[,1])]
exSTDseq <- seq[-which(names(seq) %in% unique(alignlen[,1]))]

write.fasta(STDseq, names(STDseq), "Table/03_16S_STD_seq.fasta" )
write.fasta(exSTDseq, names(exSTDseq), "Table/03_16S_exSTD_seq.fasta" )

exstdmat <- seqtab[,-which(colnames(seqtab) %in% names(STDseq))]

saveRDS(exstdmat, "Table/03_16S_excludeSTDmatrix.rds")

##############################################

## -- Assign taxonomy
s <- Sys.time()
stdmat <- seqtab[, which(colnames(seqtab) %in% names(STDseq))]
colnames(stdmat) <- sapply(STDseq[colnames(stdmat)], paste, collapse="")
stdtaxa <- assignTaxonomy(stdmat, reference, multithread=TRUE)
stdtaxa[is.na(stdtaxa)] <- "Unidentified"
print(Sys.time()-s)

stdtaxaprint <- stdtaxa
rownames(stdtaxaprint) <- names(STDseq)

stdmerege <- Taxa.mat(x=as.data.frame(seqtab[,rownames(stdtaxaprint)]), y=as.data.frame(stdtaxaprint), taxaLabel="Phylum")
saveRDS(stdmerege, "Table/03_16S_STDmatrix.rds")
