#!/bin/bash

####################################################################
# -- Vsearch v1.9.10
####################################################################

input=Table/03_ITS_exSTD_seq.fasta
out=Table
mkdir -p ${out}/code

seqtab=Table/03_ITS_excludeSTDmatrix.rds
####################################################################
for i in 0.97
do
    vsearch --cluster_fast $input --id ${i} --mothur_shared_out ${out}/04_ASV_OTU_corestab${i}.txt --centroids ${out}/04_OTUseq_${i}.fasta --msaout ${out}/04_seqAlign_${i}.txt
    
cat <<EOF > ${out}/code/OTU_aggreagate_${i}.R
## ========================================== ##
##
## Making spread sheet
##
## ========================================== ##
## -- Load library and fata tables

library(seqinr)
otu <- read.table("Table/04_ASV_OTU_corestab${i}.txt",
                  header=TRUE, row.names=2)[,-c(1:2)]
seqtab <- readRDS("Table/03_ITS_excludeSTDmatrix.rds")

## ========================================== ##

otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab),
                 dimnames=list(rownames(seqtab), colnames(otu)))

for(i in 1:ncol(otu)){
  
  if( sum(otu[,i])>1){
    memberSeq <- rownames(otu)[which(otu[,i]>0)]
    otutab[,i] <- rowSums(seqtab[,which(colnames(seqtab) %in% memberSeq) ])
  }else{
    centroidSeq <- colnames(otu)[i]
    otutab[,i] <- seqtab[, which(colnames(seqtab) == centroidSeq) ]
  }
  
}
    
saveRDS(otutab, "${out}/04_seqOTUtab_${i}.rds")

EOF

Rscript ${out}/code/OTU_aggreagate_${i}.R

done
