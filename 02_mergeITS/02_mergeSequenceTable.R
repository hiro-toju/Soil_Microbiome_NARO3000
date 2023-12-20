##############################################
##
## 2021. 5. 30
## Filtering ASV/OTU and sample
##
##############################################

# -- If you want sum up sequence read, you select "sum".
type="sum"

# -- Path to the data
input.path="runRDS"
output.path="Table"

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

# -- Loading packages
library(dada2)
library(seqinr)

## -- Loading data table
files <- list.files(input.path, full.names = TRUE)
tables=lapply(files, readRDS)

# Set working directory at containing "~stall_rmChimera.rds"
dir.create(output.path)

## ======================================= ##
## -- Combine data table
stall <- mergeSequenceTables(tables=tables, repeats=type)
seqtab.nochim <- removeBimeraDenovo(stall, method="consensus", multithread=TRUE, verbose=TRUE)

## -- Make fasta file
seq.mat <- cbind(colnames(seqtab.nochim),sprintf('X_%s', formatC(1:ncol(seqtab.nochim), width = nchar(ncol(seqtab.nochim)), flag = "0")))
write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/02_ITSmerge_clus_seq.fasta", output.path) )

seqtab <- seqtab.nochim
colnames(seqtab) <- seq.mat[,2]
saveRDS(seqtab,  sprintf("%s/02_ITSmerge_seqtab.rds", output.path))

print(sum(rowSums(stall)>2000))
## ======================================= ##


##############################################

