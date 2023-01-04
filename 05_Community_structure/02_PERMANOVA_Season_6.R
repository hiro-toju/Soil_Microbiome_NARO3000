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
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer', 'dplyr'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/02_output')

# -- Load data table
sml0 <- readRDS("Table/0301_sml.rds") 
sml2 <- read.table("./Table/SampleInfo_SamplingSeason_2.txt", header=T)

sml1 <- cbind(Sample.ID=rownames(sml0), sml0)

sml <- left_join(sml1, sml2[, c(1,4,5,6)], by="Sample.ID")
rownames(sml) <- sml$Sample.ID

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
reldf <- list("Prok"=comm16s, "Fng"=commits)

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa <- rbind(taxa16S, taxaITS)
############################################################################

per=10000
tL <- c('Family')
#tL <- c('ID', 'Genus', 'Family', 'Order')
thread=6
############################################################################

taxaL <- lapply(reldf, Taxa.mat, taxa=taxa, taxaLabel=tL)

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

res2 <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='bray',
					   parallel = thread)
	
saveRDS(res2, sprintf('%s/permanova_bray_Model1_%s_%s_season_fungi.rds', dir$rdsdir, tL, per))
sink(sprintf('%s/permanova_bray_Model1_%s_%s_season_fungi.txt', dir$tabledir, tL, per))
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

res1 <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='bray',
					   parallel = thread)
	
saveRDS(res1, sprintf('%s/permanova_bray_Model1_%s_%s_season_pro.rds', dir$rdsdir, tL, per))
sink(sprintf('%s/permanova_bray_Model1_%s_%s_season_pro.txt', dir$tabledir, tL, per))
print(res1); sink()

