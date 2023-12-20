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
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/03_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 

clr16s <- readRDS('Table/clr16s.rds')
clrits <- readRDS('Table/clrits.rds')
dflist <- list("Prok"=clr16s, "Fng"=clrits)

siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})
dflist <- lapply(dflist, function(x){ na.omit(x) })

taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################

rpca <- lapply(dflist, prcomp, scale.=TRUE, rank.=2)
saveRDS(rpca, "Table/rpca_OTU.rds")

points <- list()
points[[1]] <- rpca$Prok$x
points[[2]] <- rpca$Fng$x

names(points) <- c("Prok", "Fng")


############################################################################
# functions

col2 <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 12, name = "Set3"), palettes())

categorygg <- function(x=data, sml=sml, l=l){
    
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_manual(values=col2)+
        theme_bw(base_size=15)+
        labs(title=i, subtitle=l)
    gtmp
}

categorygg_dl <- function(x=data, sml=sml, l=l){
    
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
        theme_bw(base_size=15)+
        labs(title=i, subtitle=l)
    gtmp
}

continuousgg <- function(x=data, sml=sml, l=l){#x= points[[i]][[2]]
	
	df <- na.omit(data.frame(env=sml[rownames(x), l], x))

    gtmp <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
        theme_bw(base_size=15)+
        labs(subtitle=l)
	if(l!='CN_ratio'){
		df <- na.omit(data.frame(env=sml[rownames(x), l], x))

	    gtmp <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
	        geom_point(shape=21, alpha=0.8)+
	        scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
	        theme_bw(base_size=15)+
	        labs(subtitle=l)	
	}
    
    gtmp
}


############################################################################
## -- Visualization
############################################################################
# 'crop', 'soil_groups', 'former_crop', 'site2'

for(l in c('crop', 'soil_groups', 'former_crop', 'site2')){ #l='soil_groups'
  gtmp <- list()
    for(i in 1:2){ 
      x <- points[[i]]
      df <- cbind(sml[rownames(x), ], x)
      colnames(df)[which(colnames(df)==l)] <- 'env'
      gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
        geom_point(shape=21, alpha=0.8)+
        scale_fill_manual(values=col2)+
        theme_bw(base_size=15)+
        labs(title=i, subtitle=l) + theme(legend.position = "none")
      gtmp[[i]]     
    }
    ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
           filename=sprintf('%s/%s_clr_euclid_OTU.pdf', dir$figdir, l))
}

for(l in c('crop', 'soil_groups', 'former_crop', 'site2')){ #l='soil_groups'
  gtmp <- list()
  for(i in 1:2){ 
    x <- points[[i]]
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
      geom_point(shape=21, alpha=0.8)+
      scale_fill_manual(values=col2)+
      theme_bw(base_size=15)+
      labs(title=i, subtitle=l) + theme(legend.position = "left")
    gtmp[[i]]     
  }
  ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
         filename=sprintf('%s/%s_clr_euclid_OTU_caption.pdf', dir$figdir, l))
}

############################################################################
# 'DL'

for(l in c('DL')){ 
  gtmp <- list()
  for(i in 1:2){ 
    x <- points[[i]]
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
      geom_point(shape=21, alpha=0.8)+
      scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
      theme_bw(base_size=15)+
      labs(title=i, subtitle=l) + theme(legend.position = "none")
    gtmp[[i]]     
  }
  ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
         filename=sprintf('%s/%s_clr_euclid_OTU.pdf', dir$figdir, l))
}

for(l in c('DL')){ 
  gtmp <- list()
  for(i in 1:2){ 
    x <- points[[i]]
    df <- cbind(sml[rownames(x), ], x)
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
      geom_point(shape=21, alpha=0.8)+
      scale_fill_manual(values=c("cornflowerblue", brewer.pal(6, 'Reds')[3:6]))+
      theme_bw(base_size=15)+
      labs(title=i, subtitle=l) + theme(legend.position = "left")
    gtmp[[i]]     
  }
  ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
         filename=sprintf('%s/%s_clr_euclid_OTU_caption.pdf', dir$figdir, l))
}


############################################################################
# 'pH_dry_soil', 'available_P', 'EC_electric_conductivity', 'CNratio', 'rate_of_chemical_fertilizer_applicationN', 'rate_of_chemical_fertilizer_applicationK'

for(l in c('pH_dry_soil', 'available_P', 'EC_electric_conductivity', 'CNratio',
           'rate_of_chemical_fertilizer_applicationN', 'rate_of_chemical_fertilizer_applicationK')){ 
  gtmp <- list()
  for(i in 1:2){ 
    x <- points[[i]]
    df <- na.omit(data.frame(env=sml[rownames(x), l], x))
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
      geom_point(shape=21, alpha=0.8)+
      scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
      theme_bw(base_size=15)+
      labs(title=i) + theme(legend.position = "left")
    gtmp[[i]]     
  }
  ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
         filename=sprintf('%s/%s_clr_euclid_OTU_caption.pdf', dir$figdir, l))
}


for(l in c('pH_dry_soil', 'available_P', 'EC_electric_conductivity', 'CNratio',
           'rate_of_chemical_fertilizer_applicationN', 'rate_of_chemical_fertilizer_applicationK')){ 
  gtmp <- list()
  for(i in 1:2){ 
    x <- points[[i]]
    df <- na.omit(data.frame(env=sml[rownames(x), l], x))
    colnames(df)[which(colnames(df)==l)] <- 'env'
    gtmp[[i]] <- ggplot(df, aes(x=`PC1`, y=`PC2`, fill=env))+
      geom_point(shape=21, alpha=0.8)+
      scale_fill_gradientn(colors=brewer.pal(11, 'Spectral'), na.value = NA)+
      theme_bw(base_size=15)+
      labs(title=i) + theme(legend.position = "none")
    gtmp[[i]]     
  }
  ggsave(plot=plot_grid( plotlist=gtmp, nrow=1), w= 9, h= 4.8,
         filename=sprintf('%s/%s_clr_euclid_OTU.pdf', dir$figdir, l))
}



