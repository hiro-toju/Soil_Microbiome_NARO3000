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

randNum <- 10000
DLth <-60
octh <- 30

## -- Loading Function and Library
library(renv)
library('AnalysisHelper')
load.lib(c('igraph', "ggnetwork", 'graphlayouts',
           'RColorBrewer', 'cowplot', 'tidyr',  'dplyr'))
library(tidyverse)
library(foreach)

# -- Create directory to save
dir <- make.dir('08_Core_microbes/07_output')

# -- Load data table
read16s <- as.data.frame(readRDS('Table/rrarefy_16S.rds'))
readits <- as.data.frame(readRDS('Table/rrarefy_ITS.rds')) 

sml <- readRDS("Table/0301_sml.rds") 
siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

sml$DL[which(sml$DL=='Level 2')] <- 'Level 2-5'
sml$DL[which(sml$DL=='Level 3')] <- 'Level 2-5'
sml$DL[which(sml$DL=='Level 4')] <- 'Level 2-5'
sml$DL[which(sml$DL=='Level 5')] <- 'Level 2-5'

taxa16S <- as.data.frame(readRDS('Table/09_prok_dominatList.rds'))
taxaITS <- as.data.frame(readRDS('Table/09_fungi_dominatList.rds'))

taxa <- rbind(taxa16S, taxaITS)

############################################################################

DLrow <- rownames(sml)[which(!is.na(sml$DL))]
commonrow <- sort(intersect(intersect(rownames(read16s), rownames(readits)), DLrow), decreasing=FALSE)
dflist <- cbind(read16s[commonrow,], readits[commonrow,])
sml2 <- data.frame(sml[commonrow,], row.names = commonrow)


############################################################################
# Selecting project x crop x site2 combinations with 50 or more samples

DL.bin <- sml2$DL
DL.bin[which(DL.bin != 'Level 1')] <- 0
DL.bin[which(DL.bin == 'Level 1')] <- 1
mode(DL.bin) <- 'numeric'

df <- cbind(sml2, DL.bin)
df <- na.omit(df[, c('DL', 'DL.bin', 'experimental_purpose', 'crop', 'site2')])

df2 <- df %>% mutate(Category = paste(!!!rlang::syms(c("experimental_purpose", "crop", "site2")), sep="-"))

nCategory <- table(df2$Category)

selected <- nCategory[which(nCategory >= DLth)]
listCategory <- names(selected)

sc <- foreach(i=1:length(listCategory)) %do% {
  rownames(df2[which(df2$Category==listCategory[i]),])
}

############################################################################

matc1 <- foreach(i=1:length(listCategory)) %do% {
  dflist[sc[[i]],]
}

scol <- foreach(i=1:length(listCategory)) %do% {
  bi <- as.matrix(matc1[[i]])
  bi[which(bi > 0)] <- 1
  tbi <- t(bi)
  out <- t(subset(tbi, rowSums(tbi) >= octh))
  colnames(out)
}


matc2 <- foreach(i=1:length(listCategory)) %do% {
  dflist[sc[[i]], scol[[i]]]
}

foreach(i=1:length(listCategory)) %do% {
  saveRDS(matc2[[i]], file=sprintf('%s/Matrix_Category_%s.rds', dir$rdsdir, i))
}

sink(file=sprintf('%s/selected_list.txt', dir$tabledir))
foreach(i=1:length(listCategory), .combine = 'rbind') %do% {
  out <- c(listCategory[i], dim(matc2[[i]]))
  names(out) <- c('Category', 'N.samples', 'N.OTUs')
  out
}
sink()


diseaseLevc <- foreach(i=1:length(listCategory)) %do% {
  data.frame(disease=df2[rownames(matc2[[i]]),][,1], row.names=rownames(matc2[[i]]))
}

foreach(i=1:length(listCategory)) %do% {
  saveRDS(diseaseLevc[[i]], file=sprintf('%s/diseaseLevc_Category_%s.rds', dir$rdsdir, i))
}


###################################################################
# Randomization for each category

for.parallel = function (cores = 1)
{
  require(doParallel)
  cluster = makeCluster(cores, "FORK")
  return(registerDoParallel(cluster))
}

###################################################################

for (z in 1:length(listCategory)) {
  
  diseaseLev <- diseaseLevc[[z]]
  networkMat <- matc2[[z]]
  
  
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
  ## -- Merge randomized matrix
  
  s <- proc.time()[3]
  for.parallel(5)
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
  join[which(join$Zvalue>0),'Pvalue'] <- rowSums(orgn[which(join$Zvalue>0),'value'] < randMat[which(join$Zvalue>0),])/randNum
  join[which(join$Zvalue<0),'Pvalue'] <- rowSums(orgn[which(join$Zvalue<0),'value'] > randMat[which(join$Zvalue<0),])/randNum
  
  join[,'FDR'] <- p.adjust(join[,'Pvalue'], method='fdr')
  
  plot(join[,'FDR'], join$Zvalue)
  
  taxaD <- tail(join[order(join[,'Zvalue']), 'key'])
  taxa[as.character(taxaD),]
  
  saveRDS(join, sprintf('Table/08_07_preference_to_DL_Category_%s.rds', z))
  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
  
  wide <- spread(join[,c('level', 'key', 'Zvalue')], key=key, value=Zvalue)
  wide[is.na(wide)] <- 0
  plotorder <- hclust(dist(t(wide[,-1])))$order
  
  join$key <- factor(join$key, levels=colnames(wide)[-1][plotorder])
  
  ggDL <- ggplot(join, aes(x=level, y=key, fill=Zvalue))+
    geom_tile()+
    scale_fill_gradient2()
  ggsave(plot=ggDL, filename=sprintf('%s/Preference_to_DiseaseLevel_Category_%s.pdf', dir$figdir, z),w=5, h=15)
  
  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

}