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

#remotes::install_github("gavinsimpson/ggvegan")

set.seed(123)
setwd("../")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer', 'ggvegan'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/08_output')

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

index <- as.data.frame(readRDS('Table/0303_communityIndices_ReadCount.rds'))

sampleinfo <- data.frame(sml, index[rownames(sml), ])


############################################################################

#rpca <- lapply(dflist, prcomp, scale.=TRUE, rank.=5)
rpca <- readRDS("Table/rpca_OTU_k5.rds")

points <- list()
points[[1]] <- rpca$Prok$x
points[[2]] <- rpca$Fng$x

pcalist <- list()
pcalist[[1]] <- as.data.frame(rpca$Prok$x)
pcalist[[2]] <- as.data.frame(rpca$Fng$x)

infopro <- sampleinfo[rownames(pcalist[[1]]),c("crop", "DL", "latitude", "pH_dry_soil", "EC_electric_conductivity", "available_P", "CNratio")]
infofng <- sampleinfo[rownames(pcalist[[2]]),c("crop", "DL", "latitude", "pH_dry_soil", "EC_electric_conductivity", "available_P", "CNratio")]

dfpro <- cbind(infopro, points[[1]])
dffng <- cbind(infofng, points[[2]])

############################################################################

sink(sprintf('%s/lm_PCA_lat.txt', dir$tabledir))
summary(lm(dfpro$PC1 ~ dfpro$latitude))
summary(lm(dfpro$PC2 ~ dfpro$latitude))
summary(lm(dffng$PC1 ~ dffng$latitude))
summary(lm(dffng$PC2 ~ dffng$latitude))
sink()

############################################################################

gtmp <- list()

gtmp[[1]] <- ggplot(dfpro, aes(x=latitude, y=PC1)) +
  geom_point(shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color="red")

gtmp[[2]] <- ggplot(dfpro, aes(x=latitude, y=PC2)) +
  geom_point(shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color="red")

gtmp[[3]] <- ggplot(dffng, aes(x=latitude, y=PC1)) +
  geom_point(shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color="red")

gtmp[[4]] <- ggplot(dffng, aes(x=latitude, y=PC2)) +
  geom_point(shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color="red")


ggsave(plot=plot_grid( plotlist=gtmp, nrow=2, ncol=2), w= 6, h= 6,
       filename=sprintf('%s/Latitude_PCA.pdf', dir$figdir))

############################################################################

