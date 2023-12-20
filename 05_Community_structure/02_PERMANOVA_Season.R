############################################################################
####
#### R script for Fujita (2022)
####
#### Check community structures & PERMANOVA test
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("../")

per=1000

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer', 'dplyr'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/02_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 
clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list("Prok"=clr16s, "Fng"=clrits)

siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})
dflist <- lapply(dflist, function(x){ na.omit(x) })

taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################

tL <- c('Family')
#tL <- c('ID', 'Genus', 'Family', 'Order')
thread=6
############################################################################

taxaL <- lapply(dflist, Taxa.mat, taxa=taxa, taxaLabel=tL)

############################################################################
## Fungi

X <- taxaL[[2]]
	
smlsub <- sml[rownames(X), ]

#X <- X[1:500,]
#smlsub <- smlsub[1:500,]

expV <- c("site2", "Month", 'crop','former_crop', "soil_groups")

data.exNA <- na.omit(smlsub[,expV])

perm <- how(nperm = per)

commonrow <- intersect(rownames(data.exNA), rownames(X))
resV <- X[commonrow, which(colSums(X, na.rm=TRUE)>0)]

fml <- formula(sprintf('resV~%s', paste(expV, collapse = "+") ))

s <- proc.time()[3]
res2 <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='euclid',
					   parallel = thread)
e <- proc.time()[3]
cat(sprintf("Finished. Computation time was %.2f\n\n", e-s))

saveRDS(res2, sprintf('%s/permanova_euclid_Model1_%s_%s_season_fungi.rds', dir$rdsdir, tL, per))
sink(sprintf('%s/permanova_euclid_Model1_%s_%s_season_fungi.txt', dir$tabledir, tL, per))
print(res2); sink()


############################################################################
## Prokaryotes

X <- taxaL[[1]]
	
smlsub <- sml[rownames(X), ]
expV <- c("site2", "Month", 'crop','former_crop', "soil_groups")

data.exNA <- na.omit(smlsub[,expV])

perm <- how(nperm = per)

commonrow <- intersect(rownames(data.exNA), rownames(X))
resV <- X[commonrow, which(colSums(X, na.rm=TRUE)>0)]

fml <- formula(sprintf('resV~%s', paste(expV, collapse = "+") ))

s <- proc.time()[3]
res1 <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='euclid',
					   parallel = thread)
e <- proc.time()[3]
cat(sprintf("Finished. Computation time was %.2f\n\n", e-s))

saveRDS(res1, sprintf('%s/permanova_euclid_Model1_%s_%s_season_pro.rds', dir$rdsdir, tL, per))
sink(sprintf('%s/permanova_euclid_Model1_%s_%s_season_pro.txt', dir$tabledir, tL, per))
print(res1); sink()

