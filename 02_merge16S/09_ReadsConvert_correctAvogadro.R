############################################################################
####
#### R script for Fujita (2022)
####
#### Check standard DNA quolity
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/02_merge16S")
#### 
############################################################################

## -- Loading Function and Library
library(AnalysisHelper) 
load.lib( c('seqinr', 'ggplot2'))

# -- Load workspace and functions
stdmat <- readRDS("Table/03_16S_STDmatrix.rds")
seqtab <- readRDS('Table/07_prokaryote_seqtab.rds')
proklist <- readRDS('Table/07_prokaryote_annotation.rds')

# -- Nnumber of copies of standard DNA
std.copy.n <- c( 0.1, 0.05, 0.02, 0.01,0.005)*(6.02*10^14/1000000)/1000

############################################################################

# --  Linear regression
adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared
lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1]

############################################################################

#seqtab[seqtab<10] <- 0
seqtab <- seqtab[,colSums(seqtab)>0]

############################################################################

# -- Convert sequence reads to copy numbers
coef.summary.all <- apply(stdmat, 1, lm.coef.fun)
zero.coef <- which(!coef.summary.all > 0)
if(length(zero.coef) > 0) nc.wo.std <- c(zero.coef) # Add zero coefficient if any
coef.no.nc <- coef.summary.all[-nc.wo.std]
new.seqtab.no.nc <- as.matrix(seqtab[-nc.wo.std,])

# -- Conversion
seqtab.conv0 <- new.seqtab.no.nc/coef.no.nc

# -- Add NC samples (with raw reads)
seqtab.invalid <- seqtab[nc.wo.std,]

# -- Combine valid & invalid samples
seqtab.conv.wo.std0 <- seqtab
seqtab.conv.wo.std0[c(1:nrow(seqtab))[-nc.wo.std],] <- seqtab.conv0
seqtab.conv.wo.std0[c(1:nrow(seqtab))[nc.wo.std],] <- seqtab.invalid
seqtab.conv <- seqtab.conv.wo.std0

# -- Conversion of sample reads to calculated copy numbers
amp.factor <- 1/min(seqtab.conv[seqtab.conv != 0])
seqtab.conv.correct <- round(seqtab.conv*amp.factor) # Converted to integers

# -- Filt low quality sample (low quality means no-correlation with std.copy)
std.check <- apply(stdmat, 1, function(x){  cor( x, as.matrix(std.copy.n))  })
pdf('Figures/STDcor_totalRead.pdf'); 
plot(rowSums(stdmat), std.check); abline(h=0.9); 
plot(rowSums(stdmat), rowSums(seqtab)); abline(h= 1000); 
dev.off()
low.q.sample <- which(std.check<0.9 | rowSums(seqtab)<1000)
seqtab.conv.correct.filt <- seqtab.conv.correct[-which(rownames(seqtab.conv.correct)%in%names(low.q.sample)), ]

saveRDS(seqtab.conv.correct.filt,"Table/09_readConvert_16S.rds")
saveRDS(seqtab.conv.correct.filt,"Table/readConvert_Prok.rds")

saveRDS(seqtab[rownames(seqtab.conv.correct.filt),colnames(seqtab.conv.correct.filt)],"Table/read_Prok.rds")

rownames(proklist) <- proklist[,'ID']
proksub <- proklist[colnames(seqtab.conv.correct.filt), ]

saveRDS(proksub,"Table/09_prok_dominatList.rds")
saveRDS(proksub,"Table/taxonomyList_Prok.rds")

############################################################################

## -- Quality check
sum(seqtab)
dim(seqtab)
mean(rowSums(seqtab))
range(rowSums(seqtab))

seqtab.ex <- seqtab[-which(rownames(seqtab)%in%names(low.q.sample)), ]
length(low.q.sample)
sum(seqtab.ex)
dim(seqtab.conv.correct.filt)
mean(rowSums(seqtab.ex))
range(rowSums(seqtab.ex))

