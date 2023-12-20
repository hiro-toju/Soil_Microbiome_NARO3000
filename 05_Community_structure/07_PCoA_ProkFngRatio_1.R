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
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','SpiecEasi'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/07_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 
index <- as.data.frame(readRDS('Table/0303_communityIndices_rrarefy.rds'))


comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
reldf <- list(comm16s, commits)
#reldf <- readRDS("Table/0303_rarefiedMatrix.rds")

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################

col2 <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 12, name = "Set3"), palettes())

categorygg_woc <- function(x=data, sml=sml, l=l){
    
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp <- ggplot(df, aes(x=`Axis 1`, y=`Axis 2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_manual(values=col2)+
        theme_bw(base_size=15)+
        labs(title=i, subtitle=l) + theme(legend.position = "none")
    gtmp
}

categorygg_wocdl <- function(x=data, sml=sml, l=l){
    
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp <- ggplot(df, aes(x=`Axis 1`, y=`Axis 2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
        theme_bw(base_size=15)+
        labs(title=i, subtitle=l) + theme(legend.position = "none")
    gtmp
}

continuousgg_woc <- function(x=data, sml=sml, l=l){#x= points[[i]][[2]]
	
	df <- na.omit(data.frame(env=sml[rownames(x), l], x))

    gtmp <- ggplot(df, aes(x=`Axis.1`, y=`Axis.2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(subtitle=l) + theme(legend.position = "none")
        
	if(l!='CN_ratio'){
		df <- na.omit(data.frame(env=sml[rownames(x), l], x))

	    gtmp <- ggplot(df, aes(x=`Axis.1`, y=`Axis.2`, fill=env))+
	        geom_point(shape=21, alpha=0.8)+
	        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
	        theme_bw(base_size=15)+
	        labs(subtitle=l) + theme(legend.position = "none")
	}
  
    gtmp
}
############################################################################

continuousgg <- function(x=data, sml=sml, l=l){#x= points[[i]][[2]]
	
	df <- na.omit(data.frame(env=sml[rownames(x), l], x))

    gtmp <- ggplot(df, aes(x=`Axis.1`, y=`Axis.2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(subtitle=l)
	if(l!='CN_ratio'){
		df <- na.omit(data.frame(env=sml[rownames(x), l], x))

	    gtmp <- ggplot(df, aes(x=`Axis.1`, y=`Axis.2`, fill=env))+
	        geom_point(shape=21, alpha=0.8)+
	        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
	        theme_bw(base_size=15)+
	        labs(subtitle=l)	
	}
    
    gtmp
}
############################################################################

## -- Visualize
dimrm <- readRDS('05_Community_structure/01_output/RDS/dimensionReduction_rrarefy.rds')
pca <- lapply(dimrm, function(x){colnames(x) <- paste('Axis', 1:ncol(x)); return(x)})

points <- list(pca=pca)

############################################################################
# Without Caption



for(l in c('BbyFabundance')){ 
    glist <- c()
    for(i in names(points)){ 

        point <- lapply(points[[i]], continuousgg_woc, sml=index, l=l)
        glist[[i]] <- plot_grid( plotlist=point, nrow=1)        
    }
    ggsave(plot=plot_grid( plotlist=glist, ncol=1), w= 9, h= 4.5,
           filename=sprintf('%s/%s_bray_rrarefy_ASV.pdf', dir$figdir, l))   
}
############################################################################

############################################################################
# With Caption



for(l in c('BbyFabundance')){ 
    glist <- c()
    for(i in names(points)){ 

        point <- lapply(points[[i]], continuousgg, sml= index, l=l)
        glist[[i]] <- plot_grid( plotlist=point, nrow=1)
        
    }
    ggsave(plot=plot_grid( plotlist=glist, ncol=1), w= 11, h= 5,
           filename=sprintf('%s/%s_bray_rrarefy_ASV_caption.pdf', dir$figdir, l))
    
}
############################################################################

