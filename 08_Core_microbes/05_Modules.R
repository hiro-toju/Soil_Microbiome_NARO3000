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
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph', 'ggnetwork', 'graphlayouts', 'RColorBrewer'))
library(ggrepel)

# -- Create directory to save
dir <- make.dir('08_Core_microbes/05_output')

# -- Load data table
taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

funcCore <- readRDS('Table/0803_functionalScores.rds')
CoreL1 <- funcCore[funcCore$level=='Level 1',]; rownames(CoreL1) <- CoreL1[,'ID']

taxa [tail(CoreL1[order(CoreL1$score), 'ID'], n=10), ]
taxa [tail(CoreL1[order(CoreL1$Zvalue), 'ID'], n=10), ]

mbres <- readRDS('Table/0703_MBGraphObj.rds')
mbP <- mbres[[2]]

netInd <- readRDS('Table/0703_netInd.rds')
netIndp <- netInd[[2]]

############################################################################

all(rownames(CoreL1)==V(mbP)$name)

tbl <- data.frame(cbind(ID=V(mbP)$name, z.degree=V(mbP)$z.degree, moduleConect =V(mbP)$moduleConect))

tblcore <- CoreL1[,-c(3,4,5,6)]
summaryTab <- merge(tblcore, tbl, by="ID")
summaryTab2 <- merge(summaryTab, netIndp, by="ID")
rownames(summaryTab2) <- summaryTab2$ID

write.csv(summaryTab2, 'Table/0804_networkSummary.csv')
write.table(summaryTab2, file=sprintf('%s/0804_networkSummary.txt', dir$tabledir), sep='\t', quote=F, row.names=F)

Kingdom <- as.factor(summaryTab2$Kingdom)
Module <- as.factor(summaryTab2$module)

colvec <- palettes(18)
colvec[4] <- "gold1"
colvec[10] <- "lightpink"

############################################################################

g1 <- ggplot(summaryTab2, aes(x=betweenness, y=Zvalue)) + 
  geom_point(aes(shape=Kingdom, color=Module)) + 
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~module, ncol=4)+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/Betweenness_Zvalue_EachModule.pdf', dir$figdir), w=9, h=7)

# large size

g1 <- ggplot(summaryTab2, aes(x=betweenness, y=Zvalue)) + 
  geom_point(aes(shape=Kingdom, color=Module)) + 
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~module, ncol=4)+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/Betweenness_Zvalue_EachModule_large.pdf', dir$figdir), w=18, h=14)

############################################################################
# Olesen 2007 module analysis

# Color = Module

g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 
  geom_point(aes(shape=Kingdom, color=Module)) +
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  #geom_hline(yintercept=2.5, lty=2) + 
  #geom_vline(xintercept=0.62, lty=2) + 
  coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2))+
  facet_wrap(~module, ncol=4)+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/ModuleConnectivity_EachModule.pdf', dir$figdir), w=9, h=7)


g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 
  geom_point(aes(shape=Kingdom, color=Module)) +
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  #geom_hline(yintercept=2.5, lty=2) + 
  #geom_vline(xintercept=0.62, lty=2) + 
  coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2))+
  facet_wrap(~module, ncol=4)+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/ModuleConnectivity_EachModule_large.pdf', dir$figdir), w=18, h=14)

############################################################################
# Olesen 2007 module analysis: All

# Color = Module

g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 
  geom_point(aes(shape=Kingdom, color=Module)) +
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  #geom_hline(yintercept=2.5, lty=2) + 
  #geom_vline(xintercept=0.62, lty=2) + 
  coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2))+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/ModuleConnectivity_All.pdf', dir$figdir), w=6, h=5)

g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 
  geom_point(aes(shape=Kingdom, color=Module)) +
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  #geom_hline(yintercept=2.5, lty=2) + 
  #geom_vline(xintercept=0.62, lty=2) + 
  coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2)) +
  labs(x='Among-module connectivity', y='Within-module degree (z-standardized)')

ggsave(g1, filename=sprintf('%s/ModuleConnectivity_All_nolabel.pdf', dir$figdir), w=5, h=4)


g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 
  geom_point(aes(shape=Kingdom, color=Module)) +
  geom_tile()+ scale_shape_manual(values=c(16,17,15)) + 
  scale_color_manual(values= colvec) + 
  #geom_hline(yintercept=2.5, lty=2) + 
  #geom_vline(xintercept=0.62, lty=2) + 
  coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2))+
  geom_text_repel(aes(label =rownames(summaryTab2)), size=2, segment.size=0.25)

ggsave(g1, filename=sprintf('%s/ModuleConnectivity_All_large.pdf', dir$figdir), w=10, h=10)


