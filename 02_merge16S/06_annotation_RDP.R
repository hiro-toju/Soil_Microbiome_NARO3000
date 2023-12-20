
########################################################################
## -- DADA2 v.1.16
####################################################################
dir.create("Table/06_annotation_taxonomy", showWarnings = FALSE)

reference1 <- "Table/silva_nr99_v138.1_train_set.fa"
reference2 <- "Table/silva_species_assignment_v138.1.fa.gz"

## ===================================================== ##
## -- Load packages
library(dada2)
library(seqinr)

## -- Load data tables
seqfasta <- read.fasta("Table/04_OTUseq_0.98.fasta")
seqvec <- sapply(seqfasta, paste, collapse="")
## ===================================================== ##
# -- Assign taxonomy
taxa <- assignTaxonomy(seqvec, reference1, multithread=TRUE)
taxaAddSp <- addSpecies(taxa, reference2)

taxa.print <- taxaAddSp  # Removing sequence rownames for display only
rownames(taxa.print) <- names(seqvec)

## ===================================================== ##
# -- Compile table

# -- Compile table
if( max( sapply(strsplit(taxa.print[,"Kingdom"], "k__"), length), na.rm=TRUE ) > 1){
  taxa.print[,"Kingdom"] <- sapply(strsplit(taxa.print[,"Kingdom"], "k__"),  "[" ,2)
}

taxa.print[,"Phylum"] <- gsub("p__", "", taxa.print[,"Phylum"])
taxa.print[,"Class"] <- gsub("c__", "", taxa.print[,"Class"])
taxa.print[,"Order"] <- gsub("o__", "", taxa.print[,"Order"])
taxa.print[,"Family"] <- gsub("f__", "", taxa.print[,"Family"])
taxa.print[,"Genus"] <- gsub("g__", "", taxa.print[,"Genus"])
taxa.print[,"Species"] <- gsub("s__", "", taxa.print[,"Species"])

taxa.print[is.na(taxa.print)] <- "Unidentified"
taxa.print <- gsub("unidentified", "Unidentified", taxa.print)

write.table(cbind(ID=rownames(taxa.print) ,taxa.print), "Table/06_annotation_taxonomy/taxonomy_list.txt", row.names=FALSE, sep="\t", quote=FALSE)
saveRDS(cbind(ID=rownames(taxa.print) ,taxa.print), "Table/06_annotation_taxonomy/taxonomy_list.rds")
