
df <- read.table("Table/05_annotation_taxonomy/taxalist_NNC.tsv", header=T, row.names=1, sep="\t")
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

