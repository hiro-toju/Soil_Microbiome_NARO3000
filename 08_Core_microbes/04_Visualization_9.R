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
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph', 'ggnetwork', 'graphlayouts', 'RColorBrewer'))

# -- Create directory to save
dir <- make.dir('08_Core_microbes/04_output')

# -- Load data table
taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
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
colvec[10] <- "lightpink"

g1 <- ggplot(summaryTab2, aes(x=betweenness, y=Zvalue)) + 		geom_point(aes(shape=Kingdom, color=Module)) + geom_tile()+ scale_shape_manual(values=c(16,17,15)) + scale_color_manual(values= colvec) + geom_hline(yintercept=0, lty=2) 

ggsave(g1, filename=sprintf('%s/Betweenness_Zvalue.pdf', dir$figdir), w=4, h=3)

############################################################################
# Olesen 2007 module analysis

# Color = Module

g1 <- ggplot(summaryTab2, aes(x= as.numeric(moduleConect), y= as.numeric(z.degree))) + 		geom_point(aes(shape=Kingdom, color=Module)) + geom_tile()+ scale_shape_manual(values=c(16,17,15)) + scale_color_manual(values= colvec) + geom_hline(yintercept=2.5, lty=2) + geom_vline(xintercept=0.62, lty=2) + coord_cartesian(xlim = c(-0.02, 0.82), ylim = c(-2, 4.2))

ggsave(g1, filename=sprintf('%s/ModuleConnectivity.pdf', dir$figdir), w=4, h=3)




############################################################################


#V(mbP)$name[is.na(V(mbP)$score)] <- 0
V(mbP)$module =V(mbP)$eb

king <- summaryTab2[V(mbP)$name, 'Kingdom']
king[which(king=="Archaea")] <- 16
king[which(king=="Bacteria")] <- 17
king[which(king=="Fungi")] <- 15
mode(king) <- "numeric"

bb <- layout_as_backbone(mbP, keep=0.2)
E(mbP)$bb <- 'back'
E(mbP)$bb[bb$backbone] <- 'bone'
lay <- layout_with_mds(mbP)


#########
# Kingdom

ggnet <- ggnetwork(mbP, layout=lay)

gobj <- ggplot(ggnet, aes(x=x, y=y, xend=xend, yend=yend))+
        geom_edges(color='grey80', size=0.3)+
        #geom_edges(data=X[X$bb=='bone',], color='green3', size=1)+
        geom_nodes(aes(x=x, y=y, shape=as.factor(kingdom), color=as.factor(module)), size=4)+     
        theme_void()+
        theme(legend.box = "horizontal")+
        coord_equal()+
        scale_size(breaks=seq(-200, 400, 50), range = c(2,  8))+
        scale_shape_manual(values=c(16, 17, 15))+
        #scale_fill_manual(values=colvec)+
        scale_color_manual(values=colvec)+
        labs(fill='Module ID')+
        guides(fill = guide_legend(override.aes = list(size=4)))
        
ggsave(plot=gobj, filename=sprintf('%s/Network_Kingdom.pdf', dir$figdir), w=8, h=8)


#########
# betweenness

V(mbP)$score <- summaryTab2[V(mbP)$name, 'betweenness']
V(mbP)$score[V(mbP)$score<0] <- 0

ggnet <- ggnetwork(mbP, layout=lay)

gobj <- ggplot(ggnet, aes(x=x, y=y, xend=xend, yend=yend))+
        geom_edges(color='grey80', size=0.3)+
        #geom_edges(data=X[X$bb=='bone',], color='green3', size=1)+
        geom_nodes(aes(x=x, y=y, fill=as.factor(module), size=score), shape=21)+
        scale_fill_manual(values=colvec)+
        theme_void()+
        theme(legend.box = "horizontal")+
        coord_equal()+
        scale_size(breaks=seq(-200, 400, 50), range = c(2,  8))+
        labs(fill='Module ID')+
        guides(fill = guide_legend(override.aes = list(size=4)))

ggsave(plot=gobj, filename=sprintf('%s/Network_betweenness.pdf', dir$figdir), w=8, h=8)

##############
# eigenvector centrality

V(mbP)$score <- summaryTab2[V(mbP)$name, 'eigen']
V(mbP)$score[V(mbP)$score<0] <- 0

ggnet <- ggnetwork(mbP, layout=lay)

gobj <- ggplot(ggnet, aes(x=x, y=y, xend=xend, yend=yend))+
        geom_edges(color='grey80', size=0.3)+
        #geom_edges(data=X[X$bb=='bone',], color='green3', size=1)+
        geom_nodes(aes(x=x, y=y, fill=as.factor(module), size=score), shape=21)+
        scale_fill_manual(values=colvec)+
        theme_void()+
        theme(legend.box = "horizontal")+
        coord_equal()+
        scale_size(breaks=seq(-200, 400, 50), range = c(2,  8))+
        labs(fill='Module ID')+
        guides(fill = guide_legend(override.aes = list(size=4)))

ggsave(plot=gobj, filename=sprintf('%s/Network_eigen.pdf', dir$figdir), w=8, h=8)

##############
# Zvalue negative

V(mbP)$score <- -summaryTab2[V(mbP)$name, 'Zvalue']
V(mbP)$score[V(mbP)$score<0] <- 0

ggnet <- ggnetwork(mbP, layout=lay)

gobj <- ggplot(ggnet, aes(x=x, y=y, xend=xend, yend=yend))+
        geom_edges(color='grey80', size=0.3)+
        geom_nodes(aes(x=x, y=y, fill=as.factor(module), size=score), shape=21)+
        scale_fill_manual(values=colvec)+
        theme_void()+
        theme(legend.box = "horizontal")+
        coord_equal()+
        #scale_size(breaks=seq(-200, 400, 50), range = c(2,  8))+
        labs(fill='Module ID')+
        guides(fill = guide_legend(override.aes = list(size=4)))

ggsave(plot=gobj, filename=sprintf('%s/Network_ZvalueNega.pdf', dir$figdir), w=7, h=7)

##############
# Zvalue positive

V(mbP)$score <- summaryTab2[V(mbP)$name, 'Zvalue']
V(mbP)$score[V(mbP)$score<0] <- 0

ggnet <- ggnetwork(mbP, layout=lay)

gobj <- ggplot(ggnet, aes(x=x, y=y, xend=xend, yend=yend))+
        geom_edges(color='grey80', size=0.3)+
        geom_nodes(aes(x=x, y=y, fill=as.factor(module), size=score), shape=21)+
        scale_fill_manual(values=colvec)+
        theme_void()+
        theme(legend.box = "horizontal")+
        coord_equal()+
        #scale_size(breaks=seq(-200, 400, 50), range = c(2,  8))+
        labs(fill='Module ID')+
        guides(fill = guide_legend(override.aes = list(size=4)))

ggsave(plot=gobj, filename=sprintf('%s/Network_ZvaluePosi.pdf', dir$figdir), w=7, h=7)

############################################################################

ggDL <- ggplot(summaryTab2, aes(x=factor(module, levels=rev(1:11)), y=Zvalue, fill=as.factor(module)))+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = c(4, -4), linetype=2)+
    geom_boxplot()+
    geom_jitter()+scale_fill_manual(values= colvec)+scale_color_manual(values= colvec) + coord_flip() + 
    labs(x='Module ID', y='Specificity to disease level 1')
   
ggsave(plot=ggDL, filename=sprintf('%s/Specificity_to_DL1.pdf', dir$figdir), w=6, h=8)

ggFS <- ggplot(summaryTab2, aes(x=as.factor(module), y= score))+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = c(4, -4), linetype=2)+
    geom_boxplot()+
    geom_jitter()+
    labs(x='Module ID', y='Functional coreness')
   
ggsave(plot=ggDL, filename=sprintf('%s/Functional_coreness.pdf', dir$figdir), w=12, h=5)


ggZF <- ggplot(CoreL1, aes(x=Zvalue, y=FDR))+
        geom_vline(xintercept = c(4.0, -4), linetype=2)+
        geom_point()+
        theme_bw()+
        labs(x='Z value', y='FDR')

ggsave(plot=ggZF, filename=sprintf('%s/Zval_FDR.pdf', dir$figdir), w=7, h=5)
