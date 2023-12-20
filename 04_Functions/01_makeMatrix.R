############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### 
#### 
############################################################################

set.seed(123)
setwd("../")

# -- Create directory to save
library(AnalysisHelper)
library(seqinr)
dir <- make.dir('04_Functions/01_output')

# -- Load data table

clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list(clr16s, clrits)

taxa16S <- readRDS('Table/07_prokaryote_annotation.rds') 
seq16s <- read.fasta('Table/04_OTUseq_0.98.fasta') 
#read16s <- readRDS('Table/07_prokaryote_seqtab.rds')[rownames(dflist[[1]]),colnames(dflist[[1]])]

allits <- data.frame(readRDS('Table/07_fungi_annotation.rds'))
taxaITS <- allits[allits$ID %in% colnames(dflist[[2]]),]

############################################################################

## -- For picrust2
sub16s <- taxa16S[which(taxa16S[,1]%in% colnames(dflist[[1]])), ]
seqsub <- seq16s[rownames(sub16s)]
write.fasta(seqsub, sub16s[,1], sprintf("%s/picrutsInputFasta.fasta", dir$tabledir))

t16s <- t(dflist[[1]])
#write.table(cbind(ASV=rownames(t16s), t16s),
#            sprintf("%s/picrutsInput.txt", dir$tabledir), 
#           row.names=FALSE, quote=FALSE, sep="\t")

## -- For FAPROTAX
f16s <- taxa16S
rownames(f16s) <- f16s[,1]
f16s[,2] <- apply(f16s[,1:2], 1, paste, collapse=' ')
f16s <- gsub(' ', '_', f16s)
df <- data.frame(ASV=apply(f16s[rownames(t16s),-1], 1, paste, collapse=';'), t16s)
write.table(df,
            sprintf("%s/FAPROTAX_Input.tsv", dir$tabledir), 
            row.names=FALSE, quote=FALSE, sep="\t")


## -- For FunGUILD
#tits <- t(dflist[[2]])

taxaITSsub <- taxaITS[,-1]
initial <- c("k__",'p__','c__','o__','f__','g__','s__')
for(i in 1:ncol(taxaITSsub)){
    taxaITSsub[,i] <- paste(initial[i], taxaITSsub[,i],sep="")
}
taxainfo <- apply(taxaITSsub, 1, paste, collapse=";")
funguild <- cbind(ID=taxaITS$ID, taxonomy=taxainfo)

write.table(funguild,
            "Table/funguildInput.txt", 
            row.names=FALSE, quote=FALSE, sep="\t")

