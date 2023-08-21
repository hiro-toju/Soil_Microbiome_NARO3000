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
load.lib(c('tidyr', 'dplyr'))
library("igraph")

# -- Create directory to save
dir <- make.dir('08_Core_microbes/02_output')

# -- Load data table
taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxa <- rbind(taxa16S, taxaITS)

prokFunc <- readRDS('Table/04_06_FAPROTAX.rds')
#fngFunc <- read.table('04_FunctionAnalysis/Table/funguildInput.taxa.guilds.txt', sep='\t', header=TRUE)
fngFunc <- read.csv('Table/FungalTraits 1.2_ver_16Dec_2020.csv', header=TRUE)

mbres <- readRDS('Table/0703_MBGraphObj.rds')
mbP <- mbres[[2]]

############################################################################
length(unique(fngFunc$GENUS))

genus <- taxaITS[V(mbP)$name[grep('Fun', V(mbP)$name)],'Genus']
unique(fngFunc[fngFunc$GENUS %in% genus,])

family <- taxaITS[V(mbP)$name[grep('Fun', V(mbP)$name)],][which(genus=='Unidentified'),]
unique(fngFunc[fngFunc$Family %in% family$Family, 8 ])
fngFunc[fngFunc$Family %in% family$Family, ][fngFunc[fngFunc$Family %in% family$Family, 8 ]=="root-associated", ]

order <- taxaITS[V(mbP)$name[grep('Fun', V(mbP)$name)],][which(genus=='Unidentified'),][which(family$Family =='Unidentified'),]
unique(fngFunc[fngFunc$Order %in% order$Order, 8 ])
############################################################################


## -- Compile format of funguild
#fnglayer <- c('Plant Pathogen', 'Endophyte', 'Ectomycorrhizal', 'Fungal Parasite')
fnglayer <- c('plant_pathogen', 'arbuscular_mycorrhizal')
guildcolum <- 'primary_lifestyle'

fngmat <- matrix(0, ncol=length(fnglayer), nrow=nrow(taxaITS))
dimnames(fngmat) <- list(rownames(taxaITS), fnglayer)

for(i in fnglayer){#i= fnglayer[1]
    
    rows <- grep(i, fngFunc[, guildcolum] )
    fngsub <- fngFunc[(rows), c('GENUS', guildcolum)]
    colnames(fngsub) <- c('Genus', 'guild')
    
    joinG <- left_join(taxaITS, fngsub,'Genus' )
    
    fngmat[which(!is.na(joinG[,'guild'])), i] <- 1
    print(dim(fngsub))
}

fngFunc2 <- data.frame(ID=rownames(fngmat), fngmat)

taxaITS [rownames(fngFunc2[fngFunc2[,3]>0, ]),]

## -- Merge Fungi layer and Prokaryotes layer
prokFunc2 <- data.frame(ID=rownames(prokFunc), Nfixation=prokFunc[,grep('fixation', colnames(prokFunc))])
funclayer <- (left_join(left_join(taxa, prokFunc2, fngFunc2,by='ID'), fngFunc2, by='ID'))
funclayer[is.na(funclayer)] <- 0

saveRDS(funclayer, 'Table/0802_functionLayer_fungaltrait.rds')

############################################################################