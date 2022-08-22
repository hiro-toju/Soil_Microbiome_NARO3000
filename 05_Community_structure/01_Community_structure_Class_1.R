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
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/01_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
dflist <- list(comm16s, commits)
dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################

## -- Merge level
taxalevel <- 'Class'
n.core <- parallel::detectCores()

############################################################################

## -- Barplot
mergelist <- lapply(dflist, Taxa.mat, taxa=taxa, taxaLabel=taxalevel)

cols <- lapply(mergelist, makeColPalette, color=palettes( pal='rainbow')[1:30],
               specificName = c('Bac_Unidentified', 'Arc_Unidentified','Fng_Unidentified'), 
               specificColor = c('grey90'))

rellist <- c()
for(i in 1:2){
    
    dftmp <- mergelist[[i]]/rowSums(mergelist[[i]])
    coltmp <- cols[[i]]
    rellist[[i]] <- dftmp
    if(any(coltmp=='grey30')){
        rellist[[i]] <- cbind(dftmp[ ,names(coltmp[-which(coltmp=='grey30')]) ], 
                              'Others'=rowSums(dftmp[ ,names(coltmp[which(coltmp=='grey30')]) ]))
    }
    
}

colsub <- lapply(cols, function(x){ #x=cols[[1]]
    
    if(any(x=='grey30')) { res <- c(x[-which(x=='grey30')], 'Others'='grey30')
    }else{res <- x}
    return(res)
})

lflist <- lapply(rellist, function(x){ #x=rellist[[1]]
    
    tmp <- cbind(sml, x[rownames(x),])
    lf <- gather(tmp, key, value, -c(1:ncol(sml)))
     
    return(lf)
})

ggbase <- function(data=data, ..., base_size=10){
            ggplot(data, ...)+
            scale_x_discrete(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0))+
            theme_minimal(base_size=base_size)+
            theme(panel.border=element_rect(fill=NA, size=0.6)) }

############################################################################
# Without caption

gglist <- lapply(1:2, function(x){ #x=1
        
	data=lflist[[x]]
	dataagg <- aggregate(data$value, by=list(env=data[,'site2'], key=data$key), sum, na.rm=TRUE)

    pal=colsub[[x]]
    dataagg $key <- factor(dataagg $key, levels=rev(names(pal)))
    dataagg $env <- factor(dataagg $env, levels=sort(unique(dataagg $env)))
        
    ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color='grey30', width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
            theme(axis.text.x = element_text(angle=45, hjust=1))+
            labs(x='', y='Proportion of sequencing reads', fill=taxalevel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=T)) + theme(legend.position = "none")
        
    return(ggg)
})
    
ggsave(plot=plot_grid(plotlist=gglist, ncol=1, align='v'),
       filename=sprintf('%s/barplot_Class_rrarefy.pdf', dir$figdir), w=6, h=5)


############################################################################
# With caption

gglist <- lapply(1:2, function(x){ #x=1
        
	data=lflist[[x]]
	dataagg <- aggregate(data$value, by=list(env=data[,'site2'], key=data$key), sum, na.rm=TRUE)

    pal=colsub[[x]]
    dataagg $key <- factor(dataagg $key, levels=rev(names(pal)))
    dataagg $env <- factor(dataagg $env, levels=sort(unique(dataagg $env)))
        
    ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color='grey30', width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
            theme(axis.text.x = element_text(angle=45, hjust=1))+
            labs(x='', y='Proportion of sequencing reads', fill=taxalevel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=T))
        
    return(ggg)
})
       
ggsave(plot=plot_grid(plotlist=gglist, ncol=1, align='v'),
       filename=sprintf('%s/barplot_Class_rrarefy_caption.pdf', dir$figdir), w=8, h=16)
           
############################################################################
############################################################################
# Class

reldf <- lapply(mergelist, function(x){ na.omit(x/rowSums(x)) })
distdf <- lapply(reldf, vegdist, method='bray')
pcalist <- lapply(distdf, cmdscale)

saveRDS(pcalist, sprintf('%s/dimensionReduction_rrarefy_Class.rds', dir$rdsdir))



############################################################################
# ASV

commonrow <- intersect(rownames(dflist[[1]]),rownames(dflist[[2]]))
bf <- na.omit(cbind(dflist[[1]][commonrow, ], dflist[[2]][commonrow, ]))
bf[bf>0] <- 1

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })
distdf <- lapply(reldf, vegdist, method='bray')
pcalist <- lapply(distdf, cmdscale)
#mdslist <- lapply(reldf, metaMDS, distance='bray', parallel=parallel::detectCores())
#tsnelist <- lapply(reldf, Rtsne, check_duplicates=FALSE)
#phate <- lapply(reldf, phate)

saveRDS(pcalist, sprintf('%s/dimensionReduction_rrarefy.rds', dir$rdsdir))

#saveRDS(list(pcalist,  tsnelist, phate), sprintf('%s/dimensionReduction_rrarefy.rds', dir$rdsdir))

#pcoa <- cmdscale(vegdist(bf, method='jaccard', binary=TRUE))
#mdslist <- metaMDS(bf, distance='jaccard', parallel=parallel::detectCores())
#tsnelist <- Rtsne(bf, check_duplicates=FALSE)
#phate <- phate(bf)
#saveRDS(list(pcoa, tsnelist, phate), sprintf('%s/dimensionReduction_binary_rrarefy.rds', dir$rdsdir))

############################################################################
