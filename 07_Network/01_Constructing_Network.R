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

setwd("../")
th <- 1000

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('SpiecEasi', 'igraph'))

# -- Create directory to save
dir <- make.dir('07_Network/01_output')

# -- Load data table
read16s <- as.data.frame(readRDS('Table/07_prokaryote_seqtab.rds'))
readits <- as.data.frame(readRDS('Table/07_fungi_seqtab.rds')) 

th16s <- subset(read16s, rowSums(read16s) > th -1)
thits <- subset(readits, rowSums(readits) > th -1)

filt16s <- t(subset(t(th16s), rowSums(t(th16s)) > 0))
filtits <- t(subset(t(thits), rowSums(t(thits)) > 0))

dim(filt16s)
dim(filtits)

commonrow <- intersect(rownames(filt16s), rownames(filtits))

dflist <- list(filt16s[commonrow, colSums(filt16s)>0], filtits[commonrow, colSums(filtits)>0])

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })

sml <- readRDS("Table/0301_sml.rds") 
siteCount <- table(sml$site2)
a <- names(which(siteCount>10))
sml <- sml[sml$site2%in%a, ]

commind <- as.data.frame(readRDS('Table/0303_communityIndices_ReadCount.rds'))

############################################################################

# -- Set parameter
n.core <- 5

############################################################################

filt <- function(x, th=10, nsample=nrow(x),
                 minth=0.05, maxth=1){
    freq <- colSums(x>th)/nsample
    cat( sprintf('ASVs frequency range : minimum is %s, maximum is %s\n', 
                 round(min(freq), 2), round(max(freq), 2)))
    cat(sprintf('Discard ASVs in which appeared less than %s samples\n', round(nsample*minth)))
    y <- x[, freq>minth & freq<maxth]
    z <- y[rowSums(y)>0, ]
    cat( sprintf('Remained sample/asvs is %s/%s\n', nrow(z), ncol(z)))
    return(z) }

filtlist <- lapply(dflist, filt, minth=0.1, th=10, maxth=1)

#slrlist <- c()
#s2 <- Sys.time()
#for(rth in c(10:30)){ #rth=5
    
#    s <- Sys.time()

#    est <-  multi.spiec.easi(list(as.matrix(filtlist[[1]][commonrow,]), as.matrix(filtlist[[2]][commonrow,])), r=rth, method='slr',
#                       nlambda=50, lambda.min.ratio=1e-2,
#                       sel.criterion='bstars', pulsar.select=TRUE, 
#                       pulsar.params=list(rep.num= 100, ncores= n.core))
#    slrlist[[sprintf('%s', rth)]] <- est                   
#    e <- Sys.time()
#    print(e-s)
#}
#e2 <- Sys.time()
#print(e2-s2)

s3 <- Sys.time()
mblist <-  multi.spiec.easi(list(as.matrix(filtlist[[1]][commonrow,]), as.matrix(filtlist[[2]][commonrow,])),  
                                                        method='mb', nlambda=50, lambda.min.ratio=1e-3,
                                sel.criterion='bstars', pulsar.select=TRUE, 
                                pulsar.params=list(rep.num= 100, ncores= n.core))
e3 <- Sys.time()                               
print(e3-s3)

#saveRDS(slrlist, sprintf('%s/slrResult.rds', dir$rdsdir))                              
saveRDS(mblist, sprintf('%s/mbResult.rds', dir$rdsdir))                         

############################################################################


