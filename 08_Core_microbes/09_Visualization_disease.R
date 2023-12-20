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
load.lib(c('tidyr', 'dplyr', 'ggplot2', 'igraph', 'ggnetwork', 'graphlayouts', 'RColorBrewer'))
library(oaqc)
library(ggeffects)
library(cowplot)
library(foreach)

# -- Create directory to save
dir <- make.dir('08_Core_microbes/09_output')

# -- Load data table
table <- readRDS('Table/0806_NetworkOTUs_summary.rds')[,-c(2:5)]
x <- readRDS('Table/08_01_preference_to_DL.rds')
join <- x[order(x$key, decreasing=FALSE), ]
colnames(join)[7] <- 'Specificity'

############################################################################

mod.sample <- foreach(i=1:11) %do% {subset(table$ID, table$module==i)}
mod.join <- foreach(i=1:11) %do% {join[join$key %in% mod.sample[[i]],]}

ggDL <-  foreach(i=1:11) %do% {
  ggplot(mod.join[[i]], aes(x=level, y=key, fill=Specificity))+
    geom_tile()+
    scale_fill_gradient2()+
    theme(axis.text.y  = element_text(angle = 0, size=4), axis.title.x = element_blank(), axis.title.y = element_blank()) 
}

ggsave(plot=plot_grid(plotlist=ggDL, ncol=3, align='v'),
       filename=sprintf('%s/Specificity_to_disease_label.pdf', dir$figdir), w=15, h=20)


ggDL <-  foreach(i=1:11) %do% {
  ggplot(mod.join[[i]], aes(x=level, y=key, fill=Specificity))+
    geom_tile()+
    scale_fill_gradient2()+
    theme(axis.text.x  = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) 
}

ggsave(plot=plot_grid(plotlist=ggDL, ncol=3, align='v'),
       filename=sprintf('%s/Specificity_to_disease.pdf', dir$figdir), w=15, h=20)

