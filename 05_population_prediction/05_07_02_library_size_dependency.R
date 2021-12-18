############################################################################
####
#### R script for Fujita (2019)
####
#### Multivariate S-map predict all point
#### 2020.05 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('lmerTest', 'ggplot2', 'tidyr', 'cowplot', 'extrafont','ggsignif'))

# -- Create directory to save
dir <- make.dir('05_population_prediction/05_07_output')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/03_matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/03_04_sample_info.rds')
color <- readRDS('Table/03_color_palette.rds')

############################################################################

files <- list.files(dir$rdsdir)

############################################################################
statics.all <- pvalall <-  c()

for(i in names(dlist)[-7]){ # for test : i = names(dlist)[1]
    
    ## ======================================================= ##
    ## -- Extracting one treatment matrix
    
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i,]
    stat <- readRDS( paste(dir$rdsdir, files[grep( gsub('/','_',i),files)], sep='/'))
    
    ## ======================================================= ##
	a <- lm(smapT.rho~libsize, data=stat)
	statics.all <- rbind(statics.all, data.frame(treat=i,stat))
	
	pval <- pairwise.t.test(stat$smapT.rho,stat$libsize, p.adj='fdr', pool.sd = FALSE )
	pval.df <- as.data.frame( cbind(lib=1:6, lib2=7, value=pval$p.value[nrow(pval$p.value), ] ) )
	pval.df$sig <- 	ifelse(pval.df$value<0.05, '*', '')	
	pvalall <-  rbind(pvalall, cbind(treat=i, pval.df))
}    

statics.all$treat <- gsub('medium ', 'Medium-', statics.all$treat)
pvalall$treat <- gsub('medium ', 'Medium-', pvalall$treat)

g <- ggplot(statics.all[statics.all$tp==1, ], aes(x=as.factor(libsize), y=smapT.rho))+
    geom_jitter(aes(color=as.factor(replicate)))+
    geom_violin(alpha=0.7)+
    geom_point( stat = "summary", fun = "median", color='black')+
    geom_segment(data= pvalall, aes(x=as.factor(lib), xend=as.factor(lib2), y=(lib*0.06+1.05), yend=(lib*0.06+1.05)))+
    geom_segment(data= pvalall, aes(x=as.factor(lib), xend=as.factor(lib), y=(lib*0.06+1.02), yend=(lib*0.06+1.05)))+
    geom_segment(data= pvalall, aes(x=as.factor(lib2), xend=as.factor(lib2), y=(lib*0.06+1.02), yend=(lib*0.06+1.05)))+
    geom_text(data= pvalall, aes(x=((lib+lib2)/2) ,
                                 y=(lib*0.06+1.05), label=sig), size=6)+
    #geom_signif(comparisons=list(c('1', '7'),c('2', '7'),c('3', '7'),c('4', '7'),c('5', '7'),c('6', '7')),
    #			test='t.test', map_signif_level=F,step_increase = 0.05)+
    scale_color_manual(values=color$replicate)+
    theme_bw(base_size=20)+
    labs(x='Number of replicate communities used as reference database')+
    facet_wrap(~treat, dir='h', ncol=2)+
    theme(text=element_text(family='Arial'),
          strip.background = element_blank())

ggsave(plot=g, w=11, h=13,
       filename=sprintf('%s/libsize_depend.tiff', dir$figdir, gsub('/', '_', i)))
