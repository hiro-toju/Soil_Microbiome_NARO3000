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

randNum <- 10000

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('SpiecEasi', 'igraph', "ggnetwork", 'graphlayouts',
           'RColorBrewer', 'cowplot', 'tidyr', 'dummies', 'dplyr'))
library(tidyverse)

# -- Create directory to save
dir <- make.dir('08_Core_microbes/01_output')

# -- Load data table
read16s <- as.data.frame(readRDS('Table/rrarefy_16S.rds'))
readits <- as.data.frame(readRDS('Table/rrarefy_ITS.rds')) 

sml <- readRDS("Table/0301_sml.rds") 

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxa <- rbind(taxa16S, taxaITS)

graph <- readRDS('Table/0703_MBGraphObj.rds')

############################################################################

scale2 <- function(X){ (X-min(X,na.rm=TRUE))/(max(X,na.rm=TRUE)-min(X,na.rm=TRUE)) }


############################################################################

## -- Matrix used network analysis
DLrow <- rownames(sml)[which(!is.na(sml$DL))]
commonrow <- intersect(intersect(rownames(read16s), rownames(readits)), DLrow)
dflist <- cbind(read16s[commonrow,], readits[commonrow,])
diseaseLev <- data.frame(disease=sml[commonrow,'DL'], row.names = commonrow)

Nodes <- V(graph[[2]])$name
networkMat <- dflist[commonrow, Nodes]


## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
orgntmp <- t(Taxa.mat(t(networkMat), diseaseLev, 'disease', mean))
orgn <- gather(data.frame(level=rownames(orgntmp), orgntmp), key, value, -1 )

order <- data.frame(order=paste(orgn[,1],orgn[,2]))
orgn[,'order'] <- order

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
## -- Generating random label matrix
randlabel <- c()
for(i in 1:randNum) {
    tmp <- networkMat
    rownames(tmp) <- sample(rownames(networkMat)) 
    randlabel[[i]] <- tmp
}

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
## -- Merge randamized matrix

for.parallel = function (cores = 1)
{
    require(doParallel)
    cluster = makeCluster(cores, "FORK")
    return(registerDoParallel(cluster))
}

s <- proc.time()[3]
for.parallel(8)
Merge <- foreach::foreach(x=randlabel)%dopar%{ 
    tmp <-  t(Taxa.mat(t(x), diseaseLev, 'disease', mean))
    lf <- gather(data.frame(level=rownames(tmp), tmp), key, value, -1 )
    df <- data.frame( order=paste(lf[,1],lf[,2]), value=lf$value )
    left_join(order, df, by='order') }
e <- proc.time()[3]
print(e-s)
randMat <- matrix(NA, ncol=randNum, nrow=nrow(Merge[[1]]))
for(i in 1:length(Merge)) randMat[,i] <- Merge[[i]][,2]

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
## -- Mean and SD of random values
randStat = data.frame(order, mean=rowMeans(randMat), sd=apply(randMat, 1, sd))

## -- Z value
join <- left_join(orgn, randStat, by='order')
join[,'Zvalue'] <- (join[,'value']-join[,'mean'])/join[,'sd']

## -- P value 
join[,'Pvalue'] <- NA
join[which(join$Zvalue>0),'Pvalue'] <- rowSums(orgn[which(join$Zvalue>0),'value'] < randMat[which(join$Zvalue>0),])/1000
join[which(join$Zvalue<0),'Pvalue'] <- rowSums(orgn[which(join$Zvalue<0),'value'] > randMat[which(join$Zvalue<0),])/1000

join[,'FDR'] <- p.adjust(join[,'Pvalue'], method='fdr')

plot(join[,'FDR'], join$Zvalue)

taxaD <- tail(join[order(join[,'Zvalue']), 'key'])
taxa[as.character(taxaD),]

saveRDS(join, 'Table/08_01_preference_to_DL.rds')
## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

wide <- spread(join[,c('level', 'key', 'Zvalue')], key=key, value=Zvalue)
wide[is.na(wide)] <- 0
plotorder <- hclust(dist(t(wide[,-1])))$order

join$key <- factor(join$key, levels=colnames(wide)[-1][plotorder])

ggDL <- ggplot(join, aes(x=level, y=key, fill=Zvalue))+
        geom_tile()+
        scale_fill_gradient2()
ggsave(plot=ggDL, filename=sprintf('%s/Preference_to_DiseaseLevel.pdf', dir$figdir),w=5, h=15)

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

