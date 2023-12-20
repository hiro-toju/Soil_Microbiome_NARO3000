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
load.lib(c('SpiecEasi', 'igraph', "ggnetwork", 'graphlayouts',
           'RColorBrewer', 'cowplot', 'tidyr'))

# -- Create directory to save
dir <- make.dir('07_Network/03_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 
siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

taxa16S <- readRDS('Table/09_prok_dominatList.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/09_fungi_dominatList.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

mbres <- readRDS('07_Network/01_output/RDS/mbResult.rds')

############################################################################

calcCent <- function(x=graph){
    
    netInd <- cbind(data.frame(row.names=V(x)$name), taxa[V(x)$name, ])
    netInd$degree <- V(x)$degree <- degree(x)
    netInd$betweenness <- V(x)$betweenness <- betweenness(x, weights=NA, normalized=TRUE)
    netInd$eigen <- V(x)$eigen <- evcent(x, weights=NA)$vector
    netInd$hubscore <- V(x)$hubscore <- hub.score(x, weights=NA)$vector
    netInd$reach2 <- V(x)$readch2 <- (ego_size(x, 2)-1)/(vcount(x)-1)
    netInd$reach3 <- V(x)$readch3 <- (ego_size(x, 3)-1)/(vcount(x)-1)
    
    return(list(x, netInd))
}

calcClus <- function(x=graph){ #x=igWcent[[1]]
    
    #netInd <- cbind(data.frame(row.names=V(x)$name), taxa[V(x)$name, ])
    #info <- cluster_infomap(x, nb.trials=100, e.weights = NA)
    #netInd$infomap <- V(x)$infomap <- info$membership
    
    eb <- cluster_edge_betweenness(x, weights= NA, directed = FALSE)
    netInd$eb <- V(x)$eb <- eb$membership
    
   # louvain <- cluster_louvain(x,weights = NA)
    #netInd$louvain <- V(x)$louvain <- louvain$membership
    
    #ind <- which.max(c(modularity(info), modularity(eb), modularity(louvain)))
    #cat( sprintf('Highest modularity is %s method\n\n', c('infomap', 'edge betweenness', 'louvain')[ind] ))
    
    return(x)
}


Zdegree <- function(x){ #x=g; 
   
    module=V(x)$eb
    modNum <- unique(module)
    
    V(x)$z.degree <- rep(NA, length(V(x)))
    
    for(i in modNum){ #i=9
        
        ## ||||||||||||||||||||||||||||||||| ##
        ## -- Species within a module
        modSp <- V(x)[ module==modNum[[i]] ]
        gsub <- induced.subgraph(x, modSp)
        
        ## -- Standardized degree
        deg <- degree(gsub)
        z.degree <- (deg-mean(deg))/sd(deg)
        
        V(x)$z.degree[module==modNum[[i]]] <- z.degree
        ## ||||||||||||||||||||||||||||||||| ##
    }
    
    return(x)
    
}


degWithinModule <- function(x){
    
    module=V(x)$eb
    modNum <- unique(module)
    
    df <- data.frame(module, row.names=V(x)$name)
    adj <- as.matrix(as_adj(x))
    
    degreeToModule <- data.frame(Taxa.mat(adj, df, 'module'))

    degreeToModule2 <- degreeToModule/degree(x)
    degreeToModule3 <- degreeToModule2^2
    
    c=1-rowSums(degreeToModule3)
    
    V(x)$moduleConect <- c
    return(x)
}

############################################################################

#＃ -- Convert to igraph
mat <- mbres$est$data

sebeta <- as(symBeta(getOptBeta(mbres), mode = "ave"), 'matrix')
sebeta[as(getRefit(mbres), 'matrix')==0] <- 0
dimnames(sebeta) <- list(colnames(mat), colnames(mat))

betaP <- sebeta; betaP[betaP<0] <-0
betaN <- sebeta; betaN[betaP>0] <-0

mb <- graph_from_adjacency_matrix(sebeta, mode='undirected',
                                    weighted = TRUE, diag=FALSE)
mbP <- graph_from_adjacency_matrix(betaP, mode='undirected',
                                    weighted = TRUE, diag=FALSE)
mbN <- graph_from_adjacency_matrix(betaN, mode='undirected',
                                    weighted = TRUE, diag=FALSE)

iglist <- list(mb, mbP, mbN)
iglist <- lapply(iglist, function(x){
    V(x)$kingdom <- taxa[V(x)$name, 'Kingdom']
    return(x)
})

ig.v0 <- lapply(iglist, function(x){
    delete_vertices(x, V(x)[degree(x)==0])
})



## +++++++++++++++++++++++++++++++++ ##
## -- Centrality and community
centdf <- lapply(ig.v0, calcCent)
netInd <- lapply(centdf, function(x){ x[[2]]})
saveRDS(netInd, 'Table/0703_netInd.rds')

write.table(netInd[[2]], file=sprintf('%s/netInd.MBp.txt', dir$tabledir), row.names=F, quote=F, sep='\t')

igWcent <- lapply(centdf, function(x){ x[[1]]})

ig.withComm <- lapply(igWcent, calcClus)

igZdeg <- lapply(ig.withComm, Zdegree)
igModConnect <- lapply(igZdeg, degWithinModule)

## +++++++++++++++++++++++++++++++++ ##

pdf(sprintf('%s/degreedistribution.pdf', dir$figdir), w=4, h=4)
lapply(netInd, function(x){
    plot( sort(table(x$degree)/nrow(x)), ylab='Proportion', xlab='degree')
})
dev.off()

glist <- lapply(netInd, function(x){
    lf <- gather(x, key, value, -c(1:8))
    
    centalgg <- ggplot(lf,aes(x=Kingdom, y=value))+
                geom_boxplot()+
                geom_jitter()+
                facet_wrap(~key, scales='free')
    return(centalgg)
})
ggsave(plot= plot_grid(plotlist=glist, ncol=1), file=sprintf('%s/centarality.pdf', dir$figdir),
       w=8, h=12)	 

############################################################################
## -- Visualize


## -- layout
gnetBB <- lapply(igModConnect, function(x) { #x=ig.withComm[[1]]
    bb <- layout_as_backbone(x, keep=0.2)
    E(x)$bb <- 'back'
    E(x)$bb[bb$backbone] <- 'bone'
    ggnetwork(x, layout=bb$xy) })

gnetBBking <- lapply(gnetBB, function(X){ #X=gnetBB[[1]]
    
        gtmp <- ggplot(X, aes(x=x, y=y, xend=xend, yend=yend))+
                geom_edges(color='grey80', size=0.7)+
                #geom_edges(data=X[X$bb=='bone',], color='green3', size=1)+
                geom_nodes(aes(x=x, y=y, fill=kingdom), size=3, shape=21)+
                scale_fill_manual(values=brewer.pal(5,'Set1'), na.value=NA)+theme_void()+
                coord_equal()+
                labs(color='')
        return(gtmp)
})

gnetBBclus <- lapply(gnetBB, function(X){ #X=gnetBB[[1]]
    
    gtmp <- ggplot(X, aes(x=x, y=y, xend=xend, yend=yend))+
            geom_edges(color='grey80', size=0.7)+
            #geom_edges(data=X[X$bb=='bone',], color='green3', size=1)+
            geom_nodes(aes(x=x, y=y, fill=as.factor(eb), shape=kingdom), size=3)+
            geom_text(aes(x=x, y=y, label=as.factor(eb), shape=kingdom), size=3)+
            scale_fill_manual(values=c(palettes(), 'black'))+theme_void()+
            scale_shape_manual(values=c(21,24,22))+
            coord_equal()+
            labs(color='')
    
    return(gtmp)
})

ggsave(plot= plot_grid(plotlist=gnetBBclus, ncol=2), file=sprintf('%s/network_all_cluster.pdf', dir$figdir),
		w=14, h=15)	 
ggsave(plot= plot_grid(plotlist=gnetBBking, ncol=2), file=sprintf('%s/network_all_kingdom.pdf', dir$figdir),
       w=14, h=15)	

saveRDS(igModConnect, 'Table/0703_MBGraphObj.rds')


############################################################################
