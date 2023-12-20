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
load.lib(c('ggplot2', 'tidyr', 'vegan', "ggmap", "maps",'cowplot', 'lmerTest' ,'RColorBrewer', 'dplyr'))
library(foreach)
library(ape)
library(Biostrings)
library(ggeffects)
library(tidyverse)
library(foreach)
library(ggplot2)

source('08_Core_microbes/functions.R')

# -- Create directory to save
dir <- make.dir('08_Core_microbes/09_output')

# -- Load data table
##########################################
# Specificity scores

cat1 <- readRDS('08_Core_microbes/08_output/RDS/Specificity_Category_1.rds')
cat6 <- readRDS('08_Core_microbes/08_output/RDS/Specificity_Category_6.rds')

sub1 <- subset(cat1, cat1$Zvalue >= 6)
sub6 <- subset(cat6, cat6$Zvalue >= 6)

##########################################
# Extracting OTUs with specificity scores >= 6

mat1 <- as.matrix(readRDS('08_Core_microbes/07_output/RDS/Matrix_Category_1.rds'))
mat6 <- as.matrix(readRDS('08_Core_microbes/07_output/RDS/Matrix_Category_6.rds'))

rel1 <- mat1/rowSums(mat1)
rel6 <- mat6/rowSums(mat6)

submat1 <- rel1[, sub1$ID]
submat6 <- rel6[, sub6$ID]

##########################################
# diseaseLevc

ds1 <- readRDS('08_Core_microbes/07_output/RDS/diseaseLevc_Category_1.rds')
ds6 <- readRDS('08_Core_microbes/07_output/RDS/diseaseLevc_Category_6.rds')

ds1$disease[which(ds1$disease=="Level 1")] <- 1
ds1$disease[which(ds1$disease=="Level 2-5")] <- 0

ds6$disease[which(ds6$disease=="Level 1")] <- 1
ds6$disease[which(ds6$disease=="Level 2-5")] <- 0

mode(ds1$disease) <- "numeric"
mode(ds6$disease) <- "numeric"

df1 <- data.frame(cbind(ds1, submat1[rownames(ds1),]))
df6 <- data.frame(cbind(ds6, submat6[rownames(ds6),]))

write.table(df1, file=sprintf("%s/Data_Category_1.txt",dir$tabledir), quote=F, sep='\t')
write.table(df6, file=sprintf("%s/Data_Category_6.txt",dir$tabledir), quote=F, sep='\t')


############################################################################
# GLM

glm <- list()

glm[[1]] <- glm(disease ~ Bac_00034, family = binomial, data=df1)
glm[[2]] <- glm(disease ~ Bac_00044, family = binomial, data=df1)
glm[[3]] <- glm(disease ~ Bac_00061, family = binomial, data=df1)
glm[[4]] <- glm(disease ~ Bac_00224, family = binomial, data=df1)
glm[[5]] <- glm(disease ~ Bac_00237, family = binomial, data=df1)
glm[[6]] <- glm(disease ~ Fun_0059, family = binomial, data=df1)
glm[[7]] <- glm(disease ~ Fun_1871, family = binomial, data=df1)
glm[[8]] <- glm(disease ~ Fun_3676, family = binomial, data=df1)
glm[[9]] <- glm(disease ~ Fun_3688, family = binomial, data=df1)
glm[[10]] <- glm(disease ~ Fun_3993, family = binomial, data=df1)
glm[[11]] <- glm(disease ~ Fun_4311, family = binomial, data=df1)

glm[[12]] <- glm(disease ~ Bac_00031, family = binomial, data=df6)
glm[[13]] <- glm(disease ~ Bac_00861, family = binomial, data=df6)
glm[[14]] <- glm(disease ~ Fun_0056, family = binomial, data=df6)


gglist <- list()

gglist[[1]] <- plot(ggpredict(glm[[1]], terms="Bac_00034 [all]"), rawdata=T) + labs(title="Bac_00034", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[2]] <- plot(ggpredict(glm[[2]], terms="Bac_00044 [all]"), rawdata=T) + labs(title="Bac_00044", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[3]] <- plot(ggpredict(glm[[3]], terms="Bac_00061 [all]"), rawdata=T) + labs(title="Bac_00061", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[4]] <- plot(ggpredict(glm[[4]], terms="Bac_00224 [all]"), rawdata=T) + labs(title="Bac_00224", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[5]] <- plot(ggpredict(glm[[5]], terms="Bac_00237 [all]"), rawdata=T) + labs(title="Bac_00237", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[6]] <- plot(ggpredict(glm[[6]], terms="Fun_0059 [all]"), rawdata=T) + labs(title="Fun_0059", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[7]] <- plot(ggpredict(glm[[7]], terms="Fun_1871 [all]"), rawdata=T) + labs(title="Fun_1871", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[8]] <- plot(ggpredict(glm[[8]], terms="Fun_3676 [all]"), rawdata=T) + labs(title="Fun_3676", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[9]] <- plot(ggpredict(glm[[9]], terms="Fun_3688 [all]"), rawdata=T) + labs(title="Fun_3688", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[10]] <- plot(ggpredict(glm[[10]], terms="Fun_3993 [all]"), rawdata=T) + labs(title="Fun_3993", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[11]] <- plot(ggpredict(glm[[11]], terms="Fun_4311 [all]"), rawdata=T) + labs(title="Fun_4311", x=NULL, y=NULL) + theme(text=element_text(size=6))

gglist[[12]] <- plot(ggpredict(glm[[12]], terms="Bac_00031 [all]"), rawdata=T) + labs(title="Bac_00031", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[13]] <- plot(ggpredict(glm[[13]], terms="Bac_00861 [all]"), rawdata=T) + labs(title="Bac_00861", x=NULL, y=NULL) + theme(text=element_text(size=6))
gglist[[14]] <- plot(ggpredict(glm[[14]], terms="Fun_0056 [all]"), rawdata=T) + labs(title="Fun_0056", x=NULL, y=NULL) + theme(text=element_text(size=6))

ggsave(plot=plot_grid(plotlist=gglist, ncol=5, align='v') ,
       filename=sprintf('%s/OTUabundance_DL.pdf', dir$figdir), w=8, h=4)

sink(sprintf("%s/OTU_DL_GLMresults.txt", dir$tabledir ))
foreach(i=1:14) %do% {summary(glm[[i]])}
sink()


res <- foreach(i=1:14, .combine='rbind') %do% {
  summary(glm[[i]])$coefficients[2,]
}

rownames(res) <- c(sub1$ID, sub6$ID)

out  <- data.frame(ID=rownames(res), res, FDR=p.adjust(res[,4], method="fdr"))

write.table(out, file=sprintf('%s/OTUabundanec_DL_results.txt', dir$tabledir), quote=F, sep='\t', row.names = FALSE)
