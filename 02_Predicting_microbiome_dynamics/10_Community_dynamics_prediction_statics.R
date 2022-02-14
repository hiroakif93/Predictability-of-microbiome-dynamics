############################################################################
####
#### R script for Fujita (2019)
####
#### Finding expected early warnings signal by visualizing
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('../')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('ggplot2', 'tidyr', 'cowplot', 'scales', 'extrafont', 'RColorBrewer','lmerTest'))

# -- Create directory to save
dir <- make.dir('02_Predicting_microbiome_dynamics/Community_dynamics')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]

smltmp <- readRDS('Table/02_Error_community_prediction.rds')
div <- readRDS('Table/01_Diversity.rds')
sml <- cbind(smltmp, div[rownames(smltmp), ])

pred.comm <- readRDS("Table/02_Community_dynamics_prediction.rds")

color <- readRDS('Table/color_palette.rds')

############################################################################
## GLM and ANOVA 

stat.all.tmp <- sml
stat.all.tmp$replicate.id <- as.factor(stat.all.tmp$replicate.id)

sink(sprintf('%s/community_prediction_anova_main.txt', dir$tabledir))
a <- lmerTest::lmer(smapT_rel_tp_1 ~ as.factor(inoculum) + as.factor(resource) + shannon + bray_against_replicate + (1| replicate.id),
                    data= stat.all.tmp[stat.all.tmp $missing=='N',])
anova(a)
sink()


#########################################################
