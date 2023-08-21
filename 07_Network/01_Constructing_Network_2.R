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

setwd("/Users/toju/Dropbox/NARO_3000soil/NARO")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('SpiecEasi', 'igraph'))

# -- Create directory to save
dir <- make.dir('07_Network/01_output_2')

# -- Load data table
read16s <- as.data.frame(readRDS('Table/read_Prok.rds'))
readits <- as.data.frame(readRDS('Table/read_Fng.rds')) 
commonrow <- intersect(rownames(read16s), rownames(readits))

dflist <- list(read16s[commonrow, colSums(read16s)>0], readits[commonrow, colSums(readits)>0])

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })

sml <- readRDS("Table/0301_sml.rds") 
commind <- readRDS("Table/0303_communityIndices.rds") 


############################################################################

# -- Set parameter
n.core <- 20

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


