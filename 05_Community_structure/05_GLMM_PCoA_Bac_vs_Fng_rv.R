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

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr', 'vegan', "ggmap", "maps",'cowplot', 'lmerTest' ,'RColorBrewer', 'dplyr'))
library(igraph)
library(lme4)
library(lattice)
library(latticeExtra)
library(ggvegan)

# -- Create directory to save
dir <- make.dir('05_Community_structure/05_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

## -- Converted
sml <- readRDS("Table/0301_sml.rds") 
index <- as.data.frame(readRDS('Table/0303_communityIndices_ReadCount.rds'))

sampleinfo <- data.frame(sml, index[rownames(sml), ])

clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list("Prok"=clr16s, "Fng"=clrits)

dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})
dflist <- lapply(dflist, function(x){ na.omit(x) })

############################################################################

rpca <- lapply(dflist, prcomp, scale.=TRUE, rank.=5)
saveRDS(rpca, "Table/rpca_OTU_k5.rds")

pcalist <- list()
pcalist[[1]] <- as.data.frame(rpca$Prok$x)
pcalist[[2]] <- as.data.frame(rpca$Fng$x)

envprok <- sampleinfo[rownames(pcalist[[1]]),]
envfng <- sampleinfo[rownames(pcalist[[2]]),]

############################################################################


classV <- c()
for(i in 1:ncol(sml)){ classV <- c(classV, class(sml[,i])) }
continuous <- which(classV=='numeric')
factor <- which(classV=='character')

## -- Used for test missing sample less than 30% of no.samples
statdf <- sml[,continuous]

statdf <- na.omit(statdf[,apply(statdf, 2, function(X){ sum(is.na(X))})/nrow(statdf) < 0.3])

## -- Statics

stat <- sml[,c(colnames(sml)[continuous], 'experimental_purpose', 'Month', 'site2', 'crop')] ##
stat[,-c((ncol(stat)-3):ncol(stat))] <- apply(stat[,-c((ncol(stat)-3):ncol(stat))], 2, scale)

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

write.table(dfbac, file=sprintf('%s/PCA.Table.Bac.txt', dir$tabledir), quote=F, sep="\t", row.names=F)
write.table(dffng, file=sprintf('%s/PCA.Table.Fng.txt', dir$tabledir), quote=F, sep="\t", row.names=F)

saveRDS(dfbac, file=sprintf('%s/PCA.Table.Bac.rds', dir$rdsdir))
saveRDS(dffng, file=sprintf('%s/PCA.Table.Fng.rds', dir$rdsdir))


########################################################
# GLMM

r1 <- lme4::glmer(DL.bin~ PC1 + PC2 + PC3 + PC4 + PC5 + (1| experimental_purpose) + (1| Month) + (1| crop), family = binomial,  data= dfbac)

r2 <- lme4::glmer(DL.bin~ PC1 + PC2 + PC3 + PC4 + PC5 + (1| experimental_purpose) + (1| Month) + (1| crop), family = binomial,  data= dffng)

sink(sprintf('%s/05_GLMM_PA.results.txt', dir$tabledir))
summary(r1)
summary(r2)
sink()

########################################################
# Correlations between PC and environmental variables

#dfbac <- readRDS(sprintf('%s/PCA.Table.Bac.rds', dir$rdsdir))
#dffng <- readRDS(sprintf('%s/PCA.Table.Fng.rds', dir$rdsdir))

cormatbac <- cor(na.omit(dfbac[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'pH_dry_soil', 'EC_electric_conductivity', 'available_P', 'CNratio')]), use='pairwise.complete.obs', method='pearson')

figbac0 <- data.frame(cormatbac[1:5, 6:ncol(cormatbac)])
figbac <- figbac0[order(rownames(figbac0), decreasing=T),]

pdf(file=sprintf('%s/PC_vs_variables_Bac.pdf', dir$figdir), width=6, height=5)
levelplot(t(figbac),
          xlab = "", ylab = "",              # remove axis titles
          scales = list(x = list(rot = 45),  # change rotation for x-axis text
                        cex = 0.8), main='Bacteria' )      
dev.off()


cormatfng <- cor(na.omit(dffng[, c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'pH_dry_soil', 'EC_electric_conductivity', 'available_P', 'CNratio')]), use='pairwise.complete.obs', method='pearson')


figfng0 <- data.frame(cormatfng[1:5, 6:ncol(cormatfng)])
figfng <- figfng0[order(rownames(figfng0), decreasing=T),]

pdf(file=sprintf('%s/PC_vs_variables_Fng.pdf', dir$figdir), width=6, height=5)
levelplot(t(figfng),
          xlab = "", ylab = "",              # remove axis titles
          scales = list(x = list(rot = 45),  # change rotation for x-axis text
                        cex = 0.8), main='Fungi' )      
dev.off()


