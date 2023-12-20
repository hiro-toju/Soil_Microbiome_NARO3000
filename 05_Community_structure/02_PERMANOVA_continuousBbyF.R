############################################################################
####
#### R script for Fujita (2022)
####
#### Check community structures & PERMANOVA test
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("../")

per=1000

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/02_output')

# -- Load data table
clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list("Prok"=clr16s, "Fng"=clrits)

sml0 <- readRDS("Table/0301_sml.rds") 
siteCount <- table(sml0$site2)
a <- names(which(siteCount>10))
sml0 <- sml0[sml0$site2%in%a, ]

index <- as.data.frame(readRDS('Table/0303_communityIndices_ReadCount.rds'))
sml <- data.frame(sml0, index[rownames(sml0), ])

dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})
dflist <- lapply(dflist, function(x){ na.omit(x) })

taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)



############################################################################

tl <- tL <- c('Family')
#tL <- c('ID', 'Genus', 'Family', 'Order')
thread=4
############################################################################


	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- Merge matrix by taxonomy level
	taxaL <- lapply(dflist, Taxa.mat, taxa=taxa, taxaLabel=tl)
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- Take target variable ( NA sample < 30 %)
	cl <- c()
	for(j in 1:ncol(sml)) cl <- c(cl, class(sml[,j]))
	smlcont <- sml[,cl=="numeric"]
	nas <- apply(smlcont, 2, function(x){ sum(is.na(x))/length(x)})

	exptmp <- names(nas)[nas < 0.3]
	expV <- exptmp[-c(grep('PC', exptmp), grep('Total', exptmp), grep('rate_of_', exptmp), grep('richness', exptmp), grep('shannon', exptmp), grep('simpson', exptmp), grep('total', exptmp))]
	expV <- expV[-c(1:2)]

	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- PERMANOVA with only continuous variables
	continous <- lapply(dflist, function(X){ #X=dflist[[1]]
	
		cat( sprintf("Start PERMANOVA of %s.\n", substr(colnames(X)[3], 1, 3) ))
		s <- proc.time()[3]
		smlsub <- sml[rownames(X), ]
   
		data.exNA <- na.omit(smlsub[,c(expV, 'experimental_purpose')])
		data.exNA[,-ncol(data.exNA)] <- apply(data.exNA[,-ncol(data.exNA)], 2, scale)
		expV2 <- colnames(data.exNA)[-ncol(data.exNA)]
	
		perm <- how(nperm = per)
		setBlocks(perm) <- with(data.exNA, experimental_purpose)

		commonrow <- intersect(rownames(data.exNA), rownames(X))
		resV <- X[commonrow, which(colSums(X, na.rm=TRUE)>0)]
	
		fml <- formula(sprintf('resV~%s', paste(expV2, collapse = "+") ))
		res <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='euclid',
					   parallel = thread)
	
		e <- proc.time()[3]
		cat(sprintf("Finished. Computation time was %.2f\n\n", e-s))
		return(res)
	})
	saveRDS(continous, sprintf('%s/permanova_euclid_Model2_BbyF_%s_%s.rds', dir$rdsdir, tl, per))
	sink( sprintf('%s/permanova_euclid_Model2_BbyF_%s_%s.txt', dir$tabledir, tl, per))
	print(continous); sink()
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##



## ====================================================== ##

