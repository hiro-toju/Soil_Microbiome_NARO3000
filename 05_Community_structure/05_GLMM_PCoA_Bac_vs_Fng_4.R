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
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr', 'vegan', "ggmap", "maps",'cowplot', 'lmerTest' ,'RColorBrewer', 'dplyr'))
library(igraph)
library(lme4)
library(lattice)
library(latticeExtra)

# -- Create directory to save
dir <- make.dir('05_Community_structure/05_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

## -- Converted
sml0 <- readRDS("Table/0301_sml.rds") 
sml2 <- read.table("./Table/SampleInfo_SamplingSeason_2.txt", header=T)
sml1 <- cbind(Sample.ID=rownames(sml0), sml0)
sml <- left_join(sml1, sml2[, c(1,4,5,6)], by="Sample.ID")
rownames(sml) <- sml$Sample.ID
sml$PC2[sml$PC2>5] <- NA

abs16s <- as.data.frame(readRDS('Table/readConvert_Prok.rds'))
absits <- as.data.frame(readRDS('Table/readConvert_Fng.rds')) 
abslist <- list(abs16s, absits)

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
dflist <- list(comm16s, commits)
dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})

index <- as.data.frame(readRDS('Table/0303_communityIndices_rrarefy.rds'))

############################################################################

#commonrow <- intersect(rownames(dflist[[1]]),rownames(dflist[[2]]))
#bf <- na.omit(cbind(dflist[[1]][commonrow, ], dflist[[2]][commonrow, ]))
#bf[bf>0] <- 1

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })
distdf <- lapply(reldf, vegdist, method='bray')
pcalist0 <- lapply(distdf, cmdscale, k=5)
pcalist <- lapply(pcalist0, round, 5)

############################################################################

#pcalist <- readRDS('05_Community_structure/01_output/RDS/dimensionReduction_rrarefy.rds')

colnames(pcalist[[1]]) <- c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4', 'PCoA5')
colnames(pcalist[[2]]) <- c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4', 'PCoA5')

saveRDS(pcalist, file=sprintf('%s/PCoA_coordinates.rds', dir$rdsdir))

############################################################################

classV <- c()
for(i in 1:ncol(sml)){ classV <- c(classV, class(sml[,i])) }
continuous <- which(classV=='numeric')
factor <- which(classV=='character')

## -- Used for test missing sample less than 30% of no.samples
statdf <- sml[,continuous]

statdf <- na.omit(statdf[,apply(statdf, 2, function(X){ sum(is.na(X))})/nrow(statdf) < 0.3])

## -- Statics

stat <- sml[,c(colnames(sml)[continuous], 'experimental_purpose', 'Month')] ##
stat[,-c(ncol(stat)-1, ncol(stat))] <- apply(stat[,-c(ncol(stat)-1, ncol(stat))], 2, scale)

DL.bin <- sml$DL
DL.bin[which(DL.bin != 'Level 1')] <- 0
DL.bin[which(DL.bin == 'Level 1')] <- 1
mode(DL.bin) <- 'numeric'

df <- cbind(stat, index[rownames(stat), ], DL.bin)
df$bac.total.abundance <- log(df$bac.total.abundance)
df$fng.total.abundance <- log(df$fng.total.abundance)
df$bac.total.abundance[is.infinite(df$bac.total.abundance)] <- NA

#sink(sprintf('%s/abundance_glm_rrarefy.txt', dir$tabledir))
Sample.ID <- rownames(df)
param <- cbind(Sample.ID, df)

pcabac <- cbind(Sample.ID=rownames(pcalist[[1]]), (data.frame(pcalist[[1]])))
pcafng <- cbind(Sample.ID=rownames(pcalist[[2]]), (data.frame(pcalist[[2]])))

dfbac <- merge(pcabac, param, by='Sample.ID', all=FALSE)
dffng <- merge(pcafng, param, by='Sample.ID', all=FALSE)

write.table(dfbac, file=sprintf('%s/PCoA.Table.Bac.txt', dir$tabledir), quote=F, sep="\t", row.names=F)

write.table(dffng, file=sprintf('%s/PCoA.Table.Fng.txt', dir$tabledir), quote=F, sep="\t", row.names=F)

saveRDS(dfbac, file=sprintf('%s/PCoA.Table.Bac.rds', dir$rdsdir))
saveRDS(dffng, file=sprintf('%s/PCoA.Table.Fng.rds', dir$rdsdir))


########################################################
# GLMM


r1 <- lme4::glmer(DL.bin~ PCoA1 + PCoA2 + PCoA3 + PCoA4 + PCoA5 + (1| experimental_purpose) + (1| Month), family = binomial,  data= dfbac)

r2 <- lme4::glmer(DL.bin~ PCoA1 + PCoA2 + PCoA3 + PCoA4 + PCoA5 + (1| experimental_purpose) + (1| Month), family = binomial,  data= dffng)

sink(sprintf('%s/05_GLMM_PCoA.results.txt', dir$tabledir))
summary(r1)
summary(r2)
sink()

########################################################
# Correlations between PCoA and environmental variables

dfbac <- readRDS(sprintf('%s/PCoA.Table.Bac.rds', dir$rdsdir))
dffng <- readRDS(sprintf('%s/PCoA.Table.Fng.rds', dir$rdsdir))

cormatbac <- cor(na.omit(dfbac[, c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4', 'PCoA5', 'pH_dry_soil', 'EC_electric_conductivity', 'available_P', 'CNratio')]), use='pairwise.complete.obs', method='pearson')

color.ramp.length <- 14
negative.length <- 8
positive.length <- 6
cols <- c(colorRampPalette(c("blue", "white"))(negative.length),
          colorRampPalette(c("white", "red"))(positive.length))

figbac0 <- data.frame(cormatbac[1:5, 6:ncol(cormatbac)])
figbac <- figbac0[order(rownames(figbac0), decreasing=T),]

pdf(file=sprintf('%s/PCoA_vs_variables_Bac.pdf', dir$figdir), width=6, height=5)
levelplot(t(figbac),
          col.regions = cols,
          colorkey = list(col = cols, 
                          at = do.breaks(range(figbac), 
                                         color.ramp.length)),
          xlab = "", ylab = "",              # remove axis titles
          scales = list(x = list(rot = 45),  # change rotation for x-axis text
                        cex = 0.8), main='Bacteria' )      
dev.off()


cormatfng <- cor(na.omit(dffng[, c('PCoA1', 'PCoA2', 'PCoA3', 'PCoA4', 'PCoA5', 'pH_dry_soil', 'EC_electric_conductivity', 'available_P', 'CNratio')]), use='pairwise.complete.obs', method='pearson')

color.ramp.length <- 14
negative.length <- 8
positive.length <- 6
cols <- c(colorRampPalette(c("blue", "white"))(negative.length),
          colorRampPalette(c("white", "red"))(positive.length))

figfng0 <- data.frame(cormatfng[1:5, 6:ncol(cormatfng)])
figfng <- figfng0[order(rownames(figfng0), decreasing=T),]

pdf(file=sprintf('%s/PCoA_vs_variables_Fng.pdf', dir$figdir), width=6, height=5)
levelplot(t(figfng),
          col.regions = cols,
          colorkey = list(col = cols, 
                          at = do.breaks(range(figbac), 
                                         color.ramp.length)),
          xlab = "", ylab = "",              # remove axis titles
          scales = list(x = list(rot = 45),  # change rotation for x-axis text
                        cex = 0.8), main='Fungi' )      
dev.off()


