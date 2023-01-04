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
load.lib(c('ggplot2', 'tidyr', 'vegan', "ggmap", "maps",'cowplot', 'lmerTest' ,'RColorBrewer', 'dplyr'))
library(igraph)
library(lme4)
library(ggeffects)
library(tidyverse)
library(MASS)


# -- Create directory to save
dir <- make.dir('05_Community_structure/04_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

## -- Converted
sml0 <- readRDS("Table/0301_sml.rds") 
sml2 <- read.table("./Table/SampleInfo_SamplingSeason_2.txt", header=T)

sml1 <- cbind(Sample.ID=rownames(sml0), sml0)

sml <- left_join(sml1, sml2[, c(1,4,5,6)], by="Sample.ID")
rownames(sml) <- sml$Sample.ID

sml$PC2[sml$PC2>5] <- NA
index <- as.data.frame(readRDS('Table/0303_communityIndices_rrarefy.rds'))

############################################################################

classV <- c()
for(i in 1:ncol(sml)){ classV <- c(classV, class(sml[,i])) }
continuous <- which(classV=='numeric')
factor <- which(classV=='character')

## -- Used for test missing sample less than 30% of no.samples
statdf <- sml[,continuous]

statdf <- na.omit(statdf[,apply(statdf, 2, function(X){ sum(is.na(X))})/nrow(statdf) < 0.3])

## ============================================= ##

## -- Correlation between variables
cor.mat <- cor(statdf, index[rownames(statdf), ], use='pairwise.complete.obs',
               method='spearman')
#pairs(statdf)

lf.cor <- gather(cbind(variable=rownames(cor.mat), as.data.frame(cor.mat)),
				 key=key,value=value,-1)
				 
ggcor <- ggplot(lf.cor)+
		 geom_tile(aes(x= key, y= variable, fill=value), color='grey30')+
		 scale_fill_gradient2( low='deeppink3', limits=rep(max(abs(lf.cor$value)),2)*c(-1, 1) )+
		 coord_equal()+
         theme_minimal()+
		 theme(axis.text.x=element_text(angle=90, hjust=1))
pdf(sprintf('%s/correlation_between_variables_rrarefy.pdf', dir$figdir))
plot(ggcor)
dev.off()

## ============================================= ##
## -- Statics

stat <- sml[,c(colnames(sml)[continuous], 'experimental_purpose', 'Month')] ##
stat[,-c(ncol(stat)-1, ncol(stat))] <- apply(stat[,-c(ncol(stat)-1, ncol(stat))], 2, scale)


DL.bin <- sml$DL
DL.bin[which(DL.bin != 'Level 1')] <- 0
DL.bin[which(DL.bin == 'Level 1')] <- 1
mode(DL.bin) <- 'numeric'

df <- cbind(stat, index[rownames(stat), ], DL.bin)
df$bac.total.abundance <- log(df$bac.total.abundance)
df$fng.total.abundance <- log(df$fng.total.abundance)
df$bac.total.abundance[is.infinite(df$bac.total.abundance)] <- NA

#sink(sprintf('%s/abundance_glm_rrarefy.txt', dir$tabledir))
Sample.ID <- rownames(df)
out <- cbind(Sample.ID, df)
write.table(out, file=sprintf('%s/indexes_rrarefy.txt', dir$tabledir), quote=F, sep="\t", row.names=F)
saveRDS(df, file=sprintf('%s/indexes_rrarefy.rds', dir$rdsdir))
saveRDS(df, file='Table/indexes_rrarefy.rds')

# odd distribution found in Total_C,  rate_of_chemical_fertilizer_applicationN,P,K

########################################################
# GLMM

r1 <- lme4::glmer(DL.bin ~ scale(pH_dry_soil) + scale(EC_electric_conductivity) + scale(BbyFabundance) + (1| experimental_purpose) + (1| Month), family = binomial,  data=df, control=glmerControl(optimizer="nlminbwrap", calc.derivs = FALSE, optCtrl = list(maxfun = 10000, maxiter=100000)))

r2 <- lme4::glmer(DL.bin ~ scale(bac.shannon) + scale(fng.shannon) + scale(BbyFabundance) + (1| experimental_purpose) + (1| Month), family = binomial,  data=df)

sink(sprintf('%s/04_GLMM_variables.results_selected.txt', dir$tabledir))
summary(r1)
summary(r2)
sink()


############################################################################

#summary(lme4::glmer(DL.bin~ BbyFabundance + (1| experimental_purpose) + (1| Month), family = binomial,  data=df))

df2 <- na.omit(df[, c('BbyFabundance', 'DL.bin', 'experimental_purpose', 'Month')])
df3 <- matrix(NA, nrow=0, ncol=4)

for(i in 1:length(unique(df2$experimental_purpose))) {
	a <- subset(df2, df2$experimental_purpose ==unique(df2$experimental_purpose)[i])
	b <- nrow(a)
	c <- sum(as.numeric(a[, 'DL.bin']))
	d <- b-c
	if (d >= 10 && c >= 10) {
		df3 <- rbind(df3, a)
	} 
}

g <- ggplot(df3, aes(x= BbyFabundance, y=DL.bin)) + geom_point(shape=3) + facet_wrap(~ experimental_purpose, ncol=3)

pdf(sprintf('%s/BbyFabundance_vs_DL.bin.pdf', dir$figdir))
plot(g)
dev.off()


############################################################################


listrand <- unique(df3$experimental_purpose)

glmres <- list()
glm0a_p <- list()
gglist <- list()

for (i in 1:length(listrand)) {
	glmres[[i]] <- glm(DL.bin ~ BbyFabundance, family = binomial,  data=subset(df3, df3$experimental_purpose == listrand[i]))
	(glm0a_p[[i]] <- ggpredict(glmres[[i]], terms = "BbyFabundance"))
	gglist[[i]] <- plot(glm0a_p[[i]], rawdata=T)
}

ggsave(plot=plot_grid(plotlist=gglist, ncol=3, align='v'),
       filename=sprintf('%s/GLM_DLbin.pdf', dir$figdir), w=8, h=4)



sink(sprintf('%s/04_GLM_LD50.txt', dir$tabledir))
dose.p(glmres[[1]], p = c(0.1, 0.5, 0.9))
dose.p(glmres[[2]], p = c(0.1, 0.5, 0.9))
dose.p(glmres[[3]], p = c(0.1, 0.5, 0.9))
dose.p(glmres[[4]], p = c(0.1, 0.5, 0.9))
dose.p(glmres[[5]], p = c(0.1, 0.5, 0.9))
sink()

sink(sprintf('%s/04_GLM_DLbin.txt', dir$tabledir))
summary(glmres[[1]])
summary(glmres[[2]])
summary(glmres[[3]])
summary(glmres[[4]])
summary(glmres[[5]])
sink()

g <- ggplot(df3, aes(x= BbyFabundance, y=DL.bin)) + geom_point(shape=3) + facet_wrap(~ experimental_purpose, ncol=3)

pdf(sprintf('%s/BbyFabundance_vs_DL.bin.pdf', dir$figdir))
plot(g)
dev.off()

