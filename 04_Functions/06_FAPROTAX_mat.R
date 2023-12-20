############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("../")


df <- read.table('Table/OTU_functionList.txt')

df[,2] <- do.call(rbind, strsplit(df$V2, '_Bac'))[,1]

g <- igraph::graph_from_edgelist(as.matrix(df)[,1:2], directed=FALSE)

func <- unique(df[,1])
asv <- unique(df[,2])

mat <- matrix(0, ncol=length(func), nrow=length(asv))
dimnames(mat) <- list(asv, func)

for(i in 1:nrow(df)){
    
    f <- df[i,1]
    a <- df[i,2]
    mat[a, f] <- 1
    
}

saveRDS(mat, 'Table/04_06_FAPROTAX.rds')
############################################################################
