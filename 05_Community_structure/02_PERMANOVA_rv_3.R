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
setwd("/Users/Toju/Documents/NARO3000/Statistics")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/02_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
reldf <- list("Prok"=comm16s, "Fng"=commits)

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa <- rbind(taxa16S, taxaITS)
############################################################################

per=10000
tL <- c('Family')
#tL <- c('ID', 'Genus', 'Family', 'Order')
thread=4
############################################################################


		  
for(tl in tL){
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- Merge matrix by taxonomy level
	taxaL <- lapply(reldf, Taxa.mat, taxa=taxa, taxaLabel=tl)
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- PERMANOVA with only categorical variables
	category <- lapply(taxaL, function(X){ #X=reldf[[2]]
	
		cat( sprintf("Start PERMANOVA of %s.\n", substr(colnames(X)[3], 1, 3) ))
		s <- proc.time()[3]
		smlsub <- sml[rownames(X), ]
		expV <- c('crop','former_crop', "soil_groups", "site2", 
				  "sampling_date", 'experimental_purpose')
	
		data.exNA <- na.omit(smlsub[,expV])

		perm <- how(nperm = per)
		setBlocks(perm) <- with(data.exNA, experimental_purpose)

		commonrow <- intersect(rownames(data.exNA), rownames(X))
		resV <- X[commonrow, which(colSums(X, na.rm=TRUE)>0)]
	
		fml <- formula(sprintf('resV~%s', paste(expV[-c(length(expV))], collapse = "+") ))
		res <- adonis2(fml, data=data.exNA[commonrow,], by="margin",
					   permutations = perm, method='bray',
					   parallel = thread)
	
		e <- proc.time()[3]
		cat(sprintf("Finished. Computation time was %.2f\n\n", e-s))
		return(res)
	})
	saveRDS(category, sprintf('%s/permanova_bray_Model1_%s_%s.rds', dir$rdsdir, tl, per))
	sink( sprintf('%s/permanova_bray_Model1_%s_%s.txt', dir$tabledir, tl, per))
	print(category); sink()
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- Take target variable ( NA sample < 30 %)
	cl <- c()
	for(j in 1:ncol(sml)) cl <- c(cl, class(sml[,j]))
	smlcont <- sml[,cl=="numeric"]
	nas <- apply(smlcont, 2, function(x){ sum(is.na(x))/length(x)})

	exptmp <- names(nas)[nas < 0.3]
	expV <- exptmp[-c(grep('PC', exptmp), grep('Total', exptmp), grep('rate_of_', exptmp))]

	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
	## -- PERMANOVA with only continuous variables
	continous <- lapply(reldf, function(X){ #X=reldf[[1]]
	
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
					   permutations = perm, method='bray',
					   parallel = thread)
	
		e <- proc.time()[3]
		cat(sprintf("Finished. Computation time was %.2f\n\n", e-s))
		return(res)
	})
	saveRDS(continous, sprintf('%s/permanova_bray_Model2_%s_%s.rds', dir$rdsdir, tl, per))
	sink( sprintf('%s/permanova_bray_Model2_%s_%s.txt', dir$tabledir, tl, per))
	print(continous); sink()
	
	## ||||||||||||||||||||||||||||||||||||||||||||||||||| ##
}


## ====================================================== ##

