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
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph'))

source('08_Core_microbes/functions.R')

# -- Create directory to save
dir <- make.dir('08_Core_microbes/03_output')

# -- Load data table
taxa16S <- as.data.frame(readRDS('Table/09_prok_dominatList.rds'))
taxaITS <- as.data.frame(readRDS('Table/09_fungi_dominatList.rds'))
taxa <- rbind(taxa16S, taxaITS)

prefer <- readRDS('Table/08_01_preference_to_DL.rds')
head(prefer)
mbres <- readRDS('Table/0703_MBGraphObj.rds')

func <- readRDS('Table/0802_functionLayer_fungaltrait.rds')
rownames(func) <- func[,1]

############################################################################
lists <- c()
#for( i in 1){
	i=1
	## -- Prameters : alpha = effect of each function layer, beta is proportion of redundancy
	alpha=c(1,-i,1)
	beta=c(1,1,1)
	
	mbP <- mbres[[2]]
	colSums(func[V(mbP)$name,-c(1:8)] )

	b <- func.btw.final(g=mbP, func[V(mbP)$name,-c(1:8)], alpha, beta, gamma=1, delta=0)
	coreScore <- do.call(cbind, b)
	dimnames(coreScore) <- list(V(mbP)$name, c('score', colnames(func)[-c(1:8)]))
	
	dftmp <- data.frame(ID=rownames(coreScore), module=V(mbP)$eb, coreScore)

	colnames(prefer)[2] <- 'ID'

	summary <- left_join(prefer, dftmp, by='ID')
	CoreL1 <- summary[summary$level=='Level 1',]; rownames(CoreL1) <- CoreL1[,'ID']
	
#	lists[[i]] <- CoreL1
#}

############################################################################

saveRDS(CoreL1, 'Table/0803_functionalScores.rds')

