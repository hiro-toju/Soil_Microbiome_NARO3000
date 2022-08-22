############################################################################
####
#### R script for Fujita (2022)
####
#### Co-occurence Network
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('SpiecEasi', 'igraph', "ggnetwork", 'graphlayouts',
           'RColorBrewer', 'cowplot', 'tidyr'))

# -- Create directory to save
dir <- make.dir('07_Network/02_output')

# -- Load data table
df.16s <- readRDS('Table/readConvert_Prok.rds') ; dim(df.16s)
df.its <- readRDS('Table/readConvert_Fng.rds') ; dim(df.its)
sml <- readRDS("Table/0301_sml.rds") 

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxa <- rbind(taxa16S, taxaITS)

############################################################################

if(FALSE){
	slrres <- readRDS('07_Network/01_output/RDS/slrResult.rds')

	ebics <- c()
	for(i in 1:length(slrres)){
		
		slrtmp <- slrres[[i]]
		optind <- getOptInd(slrtmp)
		
		print(slrtmp$est$loglik[optind])
		ebics <- c(ebics , min(ebic( refit=getRefit(slrtmp), data=slrtmp$est$data, 
									 loglik=slrtmp$est$loglik[optind])))
	
	}
	plot(ebics)
	slrbest <- slrres[[which.min(ebics)]]
	saveRDS(slrbest, sprintf('%s/slrbest.rds', dir$rdsdir))
	saveRDS(slrbest, 'Table/0702_slrbest.rds')

}else{
	slrbest <- readRDS(sprintf('%s/slrbest.rds', dir$rdsdir))
}

############################################################################

calcCent <- function(x=graph){
    
    netInd <- cbind(data.frame(row.names=V(x)$name), taxa[V(x)$name, ])
    netInd$degree <- V(x)$degree <- degree(x)
    netInd$betweenness <- V(x)$betweenness <- betweenness(x, weights=NA)
    netInd$eigen <- V(x)$eigen <- evcent(x, weights=NA)$vector
    netInd$hubscore <- V(x)$hubscore <- hub.score(x, weights=NA)$vector
    netInd$reach2 <- V(x)$readch2 <- (ego_size(x, 2)-1)/(vcount(x)-1)
    netInd$reach3 <- V(x)$readch3 <- (ego_size(x, 3)-1)/(vcount(x)-1)
    
    return(list(x, netInd))
}

calcClus <- function(x=graph){ #x=igWcent[[1]]
    
    netInd <- cbind(data.frame(row.names=V(x)$name), taxa[V(x)$name, ])
    info <- cluster_infomap(x, nb.trials=100, e.weights = NA)
    netInd$infomap <- V(x)$infomap <- info$membership
    
    eb <- cluster_edge_betweenness(x, weights= NA, directed = FALSE)
    netInd$eb <- V(x)$eb <- eb$membership
    
    louvain <- cluster_louvain(x,weights = NA)
    netInd$louvain <- V(x)$louvain <- louvain$membership
    
    ind <- which.max(c(modularity(info), modularity(eb), modularity(louvain)))
    cat( sprintf('Highest modularity is %s method\n\n', c('infomap', 'edge betweenness', 'louvain')[ind] ))
    
    return(x)
}

calcClus2 <- function(x=graph, nspin=NULL){
    
    netInd <- cbind(data.frame(row.names=V(x)$name), taxa[V(x)$name, ])
    info <- cluster_infomap(x, nb.trials=100, e.weights = NA)
    netInd$infomap <- V(x)$infomap <- info$membership
    
    eb <- cluster_edge_betweenness(x, weights=NULL, directed = FALSE)
    netInd$eb <- V(x)$eb <- eb$membership
    
    louvain <- cluster_louvain(x,weights = NA)
    netInd$louvain <- V(x)$louvain <- louvain$membership
    
    if(is.null(nspin)){nspin <- max(length(unique(info$membership)),length(unique(info$eb)))  }
    spin <- cluster_spinglass(x, spins = nspin, stop.temp = 0.001, weights = NA)
    netInd$spin <- V(x)$spinglass <- spin$membership
    
    ind <- which.max(c(modularity(info), modularity(eb), modularity(louvain), modularity(spin)))
    cat( sprintf('Highest modularity is %s method\n\n', c('infomap', 'edge betweenness', 'louvain', 'spinglass')[ind] ))
    
    md.score <- c(modularity(info), modularity(eb), modularity(louvain), modularity(spin))
    names(md.score) <- c('infomap', 'edge betweenness', 'louvain', 'spinglass')
    return(list(x, netInd, md.score))
}

############################################################################

#＃ -- Convert to igraph
mat <- slrbest$est$data

icov <- Matrix::solve ( as( (slrbest$select$est$icov[[getOptInd(slrbest)]]), 'matrix') )
icov[!as(getRefit(slrbest), 'matrix')] <- 0
dimnames(icov) <- list(colnames(mat), colnames(mat))

igSLR <- graph_from_adjacency_matrix(icov, mode='undirected', 
                                     diag=FALSE, weighted=TRUE)
igP <- delete.edges(igSLR, E(igSLR)[E(igSLR)$weight < 0])
igN <- delete.edges(igSLR, E(igSLR)[E(igSLR)$weight > 0])
E(igN)$weight = abs(E(igN)$weight)

iglist <- list(igSLR, igP, igN)
iglist <- lapply(iglist, function(x){
    V(x)$kingdom <- taxa[V(x)$name, 'Kingdom']
    return(x)
})

## +++++++++++++++++++++++++++++++++ ##
## -- Centrality
centdf <- lapply(iglist, calcCent)
netInd <- lapply(centdf, function(x){ x[[2]]})
igWcent <- lapply(centdf, function(x){ x[[1]]})

ig.withComm <- lapply(igWcent, calcClus)

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

ig.v0 <- lapply(ig.withComm, function(x){
    delete_vertices(x, V(x)[degree(x)==0])
})

## -- layout
gnetBB <- lapply(ig.v0, function(x) { #x=ig.withComm[[1]]
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
            scale_fill_manual(values=palettes())+theme_void()+
            scale_shape_manual(values=c(21,24,22))+
            coord_equal()+
            labs(color='')
    
    return(gtmp)
})

ggsave(plot= plot_grid(plotlist=gnetBBclus, ncol=2), file=sprintf('%s/network_all_cluster.pdf', dir$figdir),
		w=14, h=15)	 
ggsave(plot= plot_grid(plotlist=gnetBBking, ncol=2), file=sprintf('%s/network_all_kingdom.pdf', dir$figdir),
       w=14, h=15)	

saveRDS(ig.v0, 'Table/0702_slrGraphObj.rds')
############################################################################
