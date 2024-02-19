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
dir <- make.dir('05_Community_structure/07_output')

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

infopro <- sampleinfo[rownames(pcalist[[1]]),c("crop", "DL", "pH_dry_soil", "EC_electric_conductivity", "available_P", "CNratio")]
infofng <- sampleinfo[rownames(pcalist[[2]]),c("crop", "DL", "pH_dry_soil", "EC_electric_conductivity", "available_P", "CNratio")]

envpro <- infopro[, -c(1,2)]
envfng <- infofng[, -c(1,2)]

colnames(envpro) <- c('pH', "EC", "P", "C/N")
colnames(envfng) <- c('pH', "EC", "P", "C/N")

enpro = envfit(rpca[[1]], envpro, permutations = 999, na.rm = TRUE)
enfng = envfit(rpca[[2]], envfng, permutations = 999, na.rm = TRUE)

############################################################################
col2 <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 12, name = "Set3"), palettes())

dfpro <- cbind(infopro, points[[1]])
dffng <- cbind(infofng, points[[2]])

colnames(dfpro)[which(colnames(dfpro)=="crop")] <- 'crop'
colnames(dffng)[which(colnames(dffng)=="crop")] <- 'crop'

arrow.pro <- as.data.frame(scores(enpro, display = "vectors"))
arrow.fng <- as.data.frame(scores(enfng, display = "vectors"))

############################################################################
# with labels

gtmp <- list()

gtmp[[1]] <- ggplot(dfpro) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=crop), shape=21, alpha=0.8) +
  coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=col2) +
  geom_segment(data = 50*arrow.pro,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "red")

gtmp[[2]] <- ggplot(dffng) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=crop), shape=21, alpha=0.8) +
  coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=col2) +
  geom_segment(data = 50*arrow.fng,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "red")

ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
       filename=sprintf('%s/PCA_envvec_label.pdf', dir$figdir))

############################################################################

############################################################################
# without labels

gtmp <- list()

gtmp[[1]] <- ggplot(dfpro) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=crop), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=col2) +
  geom_segment(data = 50*arrow.pro,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  theme(legend.position = "none")

gtmp[[2]] <- ggplot(dffng) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=crop), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=col2) +
  geom_segment(data = 50*arrow.fng,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  theme(legend.position = "none")
  
ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
       filename=sprintf('%s/PCA_envvec.pdf', dir$figdir))

############################################################################

gtmp <- list()

gtmp[[1]] <-autoplot(enpro, geom = 'label_repel')
gtmp[[2]] <-autoplot(enfng, geom = 'label_repel')

ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
       filename=sprintf('%s/autoplot.pdf', dir$figdir))


############################################################################
# with labels

gtmp <- list()

gtmp[[1]] <- ggplot(dfpro) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=DL), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))

gtmp[[2]] <- ggplot(dffng) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=DL), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))

ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
       filename=sprintf('%s/PCA_DL_label.pdf', dir$figdir))

############################################################################
# without labels

gtmp <- list()

gtmp[[1]] <- ggplot(dfpro) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=DL), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
  theme(legend.position = "none")

gtmp[[2]] <- ggplot(dffng) +
  geom_point(mapping=aes(x=PC1, y=PC2, fill=DL), shape=21, alpha=0.8) +
  #coord_fixed() +
  theme_bw(base_size=15) +
  scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
  theme(legend.position = "none")

ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
       filename=sprintf('%s/PCA_DL.pdf', dir$figdir))

