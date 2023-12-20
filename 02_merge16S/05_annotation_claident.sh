#!/bin/bash

########################################################################
## -- Claident ver 2019.5.10
# -- Working directory at `pwd`
# -- This script is generated at $timestamp
####################################################################
## -- Options
input=Table/04_OTUseq_1.fasta
maxpopposer=0.05
minsoratio=19
####################################################################

mkdir -p Table/05_annotation_taxonomy

# QC-auto search with the "overall_genus" sub-database.

clmakecachedb --blastdb=overall_genus --numthreads=32 $input  Table/05_annotation_taxonomy/cache.overall_genus

clidentseq --blastdb=Table/05_annotation_taxonomy/cache.overall_genus --numthreads=32 --hyperthreads=2 $input Table/05_annotation_taxonomy/Genus.clidentseq.txt


# Taonomic assignment with the LCA algorithm.
# Strict LCA assignment for the "overall_genus" search.

classigntax --taxdb=overall_genus Table/05_annotation_taxonomy/Genus.clidentseq.txt Table/05_annotation_taxonomy/Genus_strict.txt

# Relaxed LCA assignment for the "overall_genus" search.

classigntax --taxdb=overall_genus --maxpopposer=$maxpopposer --minsoratio=$minsoratio Table/05_annotation_taxonomy/Genus.clidentseq.txt Table/05_annotation_taxonomy/Genus_relax.txt


# Merging of the LCA/genus, relaxed-LCA/genus results.

clmergeassign --priority=descend Table/05_annotation_taxonomy/Genus_strict.txt Table/05_annotation_taxonomy/Genus_relax.txt Table/05_annotation_taxonomy/03.annotation.result.txt

tar cvjf Table/05_annotation_taxonomy/cache.overall_genus.tar.bz2 Table/05_annotation_taxonomy/cache.overall_genus
rm -r Table/05_annotation_taxonomy/cache.overall_genus


cat <<'KOF' > Table/05_annotation_taxonomy/claident_reformat.R

df <- read.table("Table/05_annotation_taxonomy/03.annotation.result.txt", header=T, row.names=1, sep="\t")
df2 <- df[, intersect( colnames(df), c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species") )]

for( i in 1:ncol(df2) ) {
  first.letter <- substr(colnames(df2)[i], 1, 1)
  colnames(df2)[i] <- paste( toupper(first.letter), substr(colnames(df2)[i], 2, nchar(colnames(df2)[i])), sep="")
}

mat <- as.matrix(df2)
mat[which(mat=="")] <- "Unidentified"
mat <- gsub(" ", "_", mat)

df3 <- as.data.frame(mat)
df3$identified <- "Unidentified"

for(i in 1:(ncol(df3)-1)){

   unident.row <- which(df3[,i]!="Unidentified")
   if( length(unident.row ) > 0){
      df3[unident.row,"identified"] <- as.character(df3[unident.row, i]   )	
   }

}

write.table(cbind(ID=rownames(df3), df3), "Table/05_annotation_taxonomy/taxonomy_list.txt", row.names=F, quote=F, sep="\t")

KOF

Rscript Table/05_annotation_taxonomy/claident_reformat.R


