############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### 
#### 
############################################################################

set.seed(123)
setwd("../")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr', 'vegan'))

# -- Create directory to save
dir <- make.dir('03_shaping_data/03_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

## -- Converted
abs16s <- as.data.frame(readRDS('Table/09_readConvert_16S.rds'))
absits <- as.data.frame(readRDS('Table/09_readConvert_ITS.rds')) 
abslist <- list(abs16s, absits)

#clr16s <- readRDS('Table/clr16s.rds')
#clrits <- readRDS('Table/clrits.rds')
#clr <- list(clr16s, clrits)
#saveRDS(clr, 'Table/clr_Matrix.rds')

sml <- readRDS("Table/0301_sml.rds") ; dim(sml)

################################################################

relative <- lapply(abslist, function(x){ x/rowSums(x) })
lapply(relative, dim)
## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
## -- Alpha diversity
alpha <- lapply(abslist, function(x){
				as.data.frame(cbind(richness = apply(x, 1, function(x){sum(x > 0)}),
							  		simpson  = diversity(x, index ='simpson'),
			 				  		shannon  = diversity(x, index ='shannon')
							  )	)	})

index <- cbind(bac=alpha[[1]][rownames(sml),], 
				   bac.total.abundance=rowSums(abs16s[rownames(sml),]),
				   fng=alpha[[2]][rownames(sml),], 
				   fng.total.abundance=rowSums(absits[rownames(sml),]) )
				   
Sample.ID <- rownames(sml)
				
indices <- cbind(index, BbyFabundance=log(index[,c('bac.total.abundance')])-log(index[,c('fng.total.abundance')]),
						BbyFrichness=(index[,c('bac.richness')])/(index[,c('fng.richness')]))

rownames(indices) <- Sample.ID

pdf(sprintf('%s/bacteria_fungi_ReadCount.pdf', dir$figdir), w=15, h=15)
pairs(indices); dev.off()

gg <- ggplot(indices, aes(x= bac.shannon, y= fng.shannon))+
	  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	  geom_density_2d(color='grey90', alpha=0.6)+
	  geom_point(alpha=0.5, color='orange', size=1)+
	  scale_fill_continuous(type = "viridis")+
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_continuous(expand = c(0, 0)) +
	  theme(legend.position='none',
	  		axis.title=element_text(size=15) )+
	  labs(x='Bacteria alpha diveristy', y='Fungi alpha diveristy')
	  
	  
ggsave(plot=gg, filename=sprintf('%s/bacteria_fungi_alphaDiv_ReadCount.pdf', dir$figdir), w=7, h=7)

saveRDS(indices, 'Table/0303_communityIndices_ReadCount.rds')
saveRDS(relative, 'Table/0303_RelativeMatrix_ReadCount.rds')
#saveRDS(relative, 'Table/0303_clrMatrix.rds')



write.table(cbind(Sample.ID, indices), file=sprintf('%s/0303_communityIndices_ReadCount.txt', dir$tabledir), quote=F, sep="\t", row.names=F)

############################################################################