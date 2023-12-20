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
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(ggplot2)
library(RColorBrewer)
library(grid)
library(reshape2)

## -- Create directory to save
dir <- make.dir('05_Community_structure/01_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 
siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]
dim(sml)

funcomm <- readRDS('Table/0303_RelativeMatrix_ReadCount.rds')[[2]]

relative <- funcomm[rownames(sml),]
rownames(relative) <- rownames(sml)
relative <- na.omit(relative)
relative[1:10, 1:10]

taxa <- readRDS('Table/0405_FungalTrait_all.rds')
rownames(taxa) <- taxa$ID

taxa$primary_lifestyle[which(taxa$primary_lifestyle == "")] <- 'unspecified'
taxa$primary_lifestyle[is.na(taxa$primary_lifestyle)] <- 'unspecified'

############################################################################
## -- Merge level
taxaLabel <- 'primary_lifestyle'
n.core <- parallel::detectCores()

############################################################################

## -- Barplot

comm2 <- Taxa.mat(relative, taxa=taxa, taxaLabel= taxaLabel)
comm3 <- comm2[, order(colSums(comm2), decreasing=T)]

col1 <- brewer.pal(7, "Set1")
col2 <- brewer.pal(7, "Set2")
col3 <- brewer.pal(9, "Set3")
#col4 <- brewer.pal(8, "Pastel1")

colvec <- c(col1, col2, col3)

cols <- makeColPalette(comm3, color= colvec,
               specificName = c('unspecified'), specificColor = c('grey90'))

cols[5] <- "lightpink" 
cols[6] <- "bisque1" 
cols[12] <- "yellow" 

    dftmp <- comm3/rowSums(comm3)
    coltmp <- cols
    rellist <- dftmp
    if(any(coltmp=='grey30')){
        rellist <- cbind(dftmp[ ,names(coltmp[-which(coltmp=='grey30')]) ], 
                              'Others'=rowSums(dftmp[ ,names(coltmp[which(coltmp=='grey30')]) ]))
    }
    


colsub <- lapply(cols, function(x){ #x=cols[[1]]
    
    if(any(x=='grey30')) { res <- c(x[-which(x=='grey30')], 'Others'='grey30')
    }else{res <- x}
    return(res)
})

funca <- function(x){ #x=rellist[[1]]
    tmp <- cbind(sml[rownames(comm3),], x[rownames(x),])
    lf <- gather(tmp, key, value, -c(1:ncol(sml)))
    return(lf)
}

lflist <- funca(rellist)



ggbase <- function(data=data, ..., base_size=10){
            ggplot(data, ...)+
            scale_x_discrete(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0))+
            theme_minimal(base_size=base_size)+
            theme(panel.border=element_rect(fill=NA, size=0.6)) }


       
	data=lflist
	dataagg <- aggregate(data$value, by=list(env=data[,'site2'], key=data$key), sum, na.rm=TRUE)

    pal=colsub
    dataagg $key <- factor(dataagg $key, levels=rev(names(pal)))
    dataagg $env <- factor(dataagg $env, levels=sort(unique(dataagg $env)))

############################################################################
# Without caption
        
    ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color='grey30', width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
            theme(axis.text.x = element_text(angle=45, hjust=1))+
            labs(x='', y='Proportion of sequencing reads', fill= taxaLabel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=T))+ theme(legend.position = "none")
        
plot(ggg)
   
ggsave(plot=plot_grid(ggg, ncol=1, align='v'),
       filename=sprintf('%s/barplot_FungalTraits.pdf', dir$figdir), w=6, h=3)
       
############################################################################
# With caption
        
    ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color='grey30', width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
            theme(axis.text.x = element_text(angle=45, hjust=1))+
            labs(x='', y='Proportion of sequencing reads', fill= taxaLabel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=F))
        
plot(ggg)
   
ggsave(plot=plot_grid(ggg, ncol=1, align='v'),
       filename=sprintf('%s/barplot_FungalTraits_caption.pdf', dir$figdir), w=10, h=6)
       




