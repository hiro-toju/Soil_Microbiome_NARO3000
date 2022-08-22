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

## ASV level ##

set.seed(123)
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")

#quartzFonts(ARIAL = c("ArialMT", "Arial-ItalicMT", "Arial-BoldMT", "Arial-BoldItalicMT"))
#par(family="ARIAL")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr', 'vegan', "ggmap", "maps",'cowplot', 'lmerTest' ,'RColorBrewer'))
library(latticeExtra)
library(protrix)

# -- Create directory to save
dir <- make.dir('05_Community_structure/06_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')

commonrow <- intersect(rownames(comm16s), rownames(commits))

index <- as.data.frame(readRDS('05_Community_structure/04_output/RDS/indexes_rrarefy.rds'))

sml <- readRDS("Table/0301_sml.rds") # Reseaerch sites with less than 10 samples were not included in "Table/0301_sml.rds"

common.commsml <- intersect(commonrow, rownames(sml))

sel16s0 <- comm16s[common.commsml, ]
selits0 <- commits[common.commsml, ]
sel16s <- sel16s0[, colSums(sel16s0) > 0]
selits <- selits0[, colSums(selits0) > 0]
selsml <- sml[common.commsml,]
selindex <- index[common.commsml,]

# cbind(rownames(sel16s), rownames(selits), rownames(selsml), rownames(selindex))

############################################################################
taxalevel <- 'Family'

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxcomm16s <-  Taxa.mat(sel16s, taxa= taxa16S, taxaLabel=taxalevel)
taxcommits <-  Taxa.mat(selits, taxa= taxaITS, taxaLabel=taxalevel)

############################################################################
# alpha diversity

shannon.16s <- diversity(sel16s, index ='shannon')
shannon.its <- diversity(selits, index ='shannon')


alpha <- data.frame(cbind(tapply(shannon.16s, selsml$site2, mean), tapply(shannon.16s, selsml$site2, sd), tapply(shannon.its, selsml$site2, mean), tapply(shannon.its, selsml$site2, sd)))

colnames(alpha) <- c('Shannon_prk', 'Shannon_prk_SD', 'Shannon_fng', 'Shannon_fng_SD')

saveRDS(alpha, sprintf('%s/Diversity_alpha_ASV.rds', dir$rdsdir))

write.table(cbind(site2=rownames(alpha), alpha), file=sprintf('%s/Diversity_alpha_ASV.txt', dir$tabledir), sep='\t', quote=F, row.names=F)




# beta diversity within sites
# sel16s

sub16s <- list()

for (i in 1:length(unique(selsml $site2))) {
	a <- subset(sel16s, selsml $site2==unique(selsml $site2)[i])
	sub16s[[i]] <- a[, colSums(a) > 0]
}

beta16s <- lapply(sub16s, vegdist, method="bray")

hist(unlist(beta16s))


subits <- list()

for (i in 1:length(unique(selsml $site2))) {
	a <- subset(selits, selsml $site2==unique(selsml $site2)[i])
	subits[[i]] <- a[, colSums(a) > 0]
}

betaits <- lapply(subits, vegdist, method="bray")

hist(unlist(betaits))


beta <- matrix(NA, nrow=length(unique(selsml $site2)), ncol=4)

for (i in 1:length(unique(selsml $site2))) {
	prk <- unlist(beta16s[[i]])
	Beta_Mean_prk <- mean(prk)
	Beta_SD_prk <- sd(prk)
	fng <- unlist(betaits[[i]])
	Beta_Mean_fng <- mean(fng)
	Beta_SD_fng <- sd(fng)
	beta[i,] <- c(Beta_Mean_prk, Beta_SD_prk, Beta_Mean_fng, Beta_SD_fng)
}

rownames(beta) <- unique(selsml $site2)

colnames(beta) <- c('Beta_Mean_prk', 'Beta_SD_prk', 'Beta_Mean_fng', 'BetaSD_fng')

saveRDS(beta, sprintf('%s/Diversity_beta_withinsites_ASV.rds', dir$rdsdir))

write.table(cbind(site2=rownames(beta), beta), file=sprintf('%s/Diversity_beta_withinsites_ASV.txt', dir$tabledir), sep='\t', quote=F, row.names=F)






# beta diversity within sites

commsite16s <- matrix(NA, nrow=length(unique(selsml $site2)), ncol=ncol(sel16s))
commsiteits <- matrix(NA, nrow=length(unique(selsml $site2)), ncol=ncol(selits))


for (i in 1:ncol(sel16s)) {
	commsite16s[,i] <- tapply(sel16s[,i], selsml$site2, mean)
	}

for (i in 1:ncol(selits)) {
	commsiteits[,i] <- tapply(selits[,i], selsml$site2, mean)
	}

bray_meta_prk <- vegdist(commsite16s, method="bray")
bray_meta_fng <- vegdist(commsiteits, method="bray")

beta_meta <- data.frame(cbind(Beta_meta_prk=unlist(bray_meta_prk), Beta_meta_fng=unlist(bray_meta_fng)))


saveRDS(beta_meta, sprintf('%s/Diversity_beta_amongsites_ASV.rds', dir$rdsdir))

write.table(beta_meta, file=sprintf('%s/Diversity_beta_amongsites_ASV.txt', dir$tabledir), sep='\t', quote=F, row.names=F)

Taxonomy <- c(rep("Prokaryotes", times=nrow(beta_meta)), rep("Fungi", times=nrow(beta_meta)))

Beta_btw_sites <- c(beta_meta[, 1], beta_meta[, 2])

df1 <- data.frame(cbind(Taxonomy, Beta_btw_sites))

t.test(x=beta_meta[, 1],y=beta_meta[, 2],var.equal=F,paired=F)

############################################################################
# Visualization

g <- list()

g[[1]] <- ggplot(as.data.frame(alpha), aes(x= Shannon_prk, y= Shannon_fng)) + 
	geom_errorbarh(aes(xmin= Shannon_prk - Shannon_prk_SD, xmax=Shannon_prk + Shannon_prk_SD), color="grey60") +
	geom_errorbar(aes(ymin= Shannon_fng - Shannon_fng_SD, ymax= Shannon_fng + Shannon_fng_SD), color="grey60") +
	geom_point(shape=16, size=2, col='deeppink3') + 
	labs(x='Shannon diversity (prokaryotes)', y='Shannon diversity (fungi)') + 
	geom_abline(intercept=0, slope=1, lty=3) +
	coord_cartesian(xlim = c(2, 5.5), ylim = c(2, 5.5))
	

g[[2]] <- ggplot(as.data.frame(beta), aes(x= Beta_Mean_prk, y= Beta_Mean_fng)) + 
	geom_errorbarh(aes(xmin= Beta_Mean_prk - Beta_SD_prk, xmax= Beta_Mean_prk + Beta_SD_prk), color="grey60") +
	geom_errorbar(aes(ymin= Beta_Mean_fng - BetaSD_fng, ymax= Beta_Mean_fng + BetaSD_fng), color="grey60") + 
	geom_point(shape=16, size=2, col='deepskyblue3') + 
	labs(x='Mean beta diversity within research site (prokaryotes)', y='Mean beta diversity within research site (fungi)') + 
	geom_abline(intercept=0, slope=1, lty=3) 


g[[3]] <- ggplot(df1) + geom_boxplot(aes(x=reorder(Taxonomy, as.numeric(Beta_btw_sites)), y=as.numeric(Beta_btw_sites), fill=Taxonomy)) + 
	labs(x='', y='Beta diversity between research sites') 
	

#g[[3]] <- ggplot(as.data.frame(beta_meta), aes(x= Beta_meta_prk, y= Beta_meta_fng)) + 
#	geom_point(shape=16, size=2, col='goldenrod1') + 
#	labs(x='Beta diversity between research sites (prokaryotes)', y='Beta diversity between research sites (fungi)') + 
#	geom_abline(intercept=0, slope=1, lty=3) 

ggsave(plot=plot_grid( plotlist=g, ncol=3), w= 12, h= 4.2,
           filename=sprintf('%s/Diversity_prk_vs_fng_ASV.pdf', dir$figdir))  



