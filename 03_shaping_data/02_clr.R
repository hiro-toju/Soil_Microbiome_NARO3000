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

th <- 1000

## -- Loading Function and Library
library('AnalysisHelper')
#install.packages("compositions")
#library(compositions)
library(ALDEx2)
library(foreach)
library(doParallel)
library(vegan)

#cores <- detectCores()-1
#cl <- makeCluster(cores)
#registerDoParallel(cl)


# -- Create directory to save
dir <- make.dir('03_shaping_data/02_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table
read16s <- as.data.frame(readRDS('Table/07_prokaryote_seqtab.rds'))
readits <- as.data.frame(readRDS('Table/07_fungi_seqtab.rds')) 

th16s <- subset(read16s, rowSums(read16s) >= th)
thits <- subset(readits, rowSums(readits) >= th)

filt16s <- t(subset(t(th16s), rowSums(t(th16s)) > 0))
filtits <- t(subset(t(thits), rowSums(t(thits)) > 0))

dim(filt16s)
dim(filtits)

saveRDS(filt16s, file="Table/filt_16s.rds")
saveRDS(filtits, file="Table/filt_its.rds")

min(rowSums(filt16s))
min(rowSums(filtits))
min(colSums(filt16s))
min(colSums(filtits))

rrarefied16s <- rrarefy(filt16s, th)
rarefy16s <- t(subset(t(rrarefied16s), rowSums(t(rrarefied16s)) > 0))
saveRDS(rarefy16s, file="Table/rrarefy_16S.rds")

rrarefiedits <- rrarefy(filtits, th)
rarefyits <- t(subset(t(rrarefiedits), rowSums(t(rrarefiedits)) > 0))
saveRDS(rarefyits, file="Table/rrarefy_ITS.rds")

################################################################
## clr: https://doi.org/10.3389/fmicb.2017.02224 

nrep <- 10

aldexits <- aldex.clr(t(filtits), rownames(filtits), mc.samples= nrep, denom="all")
x <- foreach(i=1:nrep) %do% {t(getMonteCarloSample(aldexits, i))}
clrits <- (x[[1]] + x[[2]] + x[[3]] + x[[4]] + x[[5]] + x[[6]] + x[[7]] + x[[8]] + x[[9]] + x[[10]])/nrep

saveRDS(clrits, file="Table/clrits.rds")

dim(clrits)

################################################################

aldex16s <- aldex.clr(t(filt16s), rownames(filt16s), mc.samples= nrep, denom="all")
x <- foreach(i=1:nrep) %do% {t(getMonteCarloSample(aldex16s, i))}
clr16s <- (x[[1]] + x[[2]] + x[[3]] + x[[4]] + x[[5]] + x[[6]] + x[[7]] + x[[8]] + x[[9]] + x[[10]])/nrep

saveRDS(clr16s, file="Table/clr16s.rds")

dim(clr16s)

