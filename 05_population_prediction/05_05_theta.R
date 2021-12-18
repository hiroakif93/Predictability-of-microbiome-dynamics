############################################################################
####
#### R script for Fujita (2019)
####
#### Finding expected early warnings signal by visualizing
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
source("functions/BidirectSimplex_v0.7.4.R")
source("functions/BidirectBlocklnlp_v0.7.4.R")

load.lib( c('rEDM', 'ggplot2', 'tidyr', 'cowplot'))

# -- Create directory to save
dir <- make.dir('05_population_prediction/05_05_output')

# -- Load data table
dlist  <- readRDS("Table/03_matrixList.rds")
libmat=dlist[[7]]
sml <- readRDS('Table/03_03_sample_info.rds')

# -- Decided regularized parameter range
E.range=1:20
theta_test <- c(0, 0.001, 0.01,0.05,0.1,0.2,0.5,1,2,4,8)
tp=1

for.parallel(8)

############################################################################

predlistT <- predlist0 <- statall <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ============================================= ##
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ]

    tmp <- lapply(1:8, function(x){
        y <- apply(ts[libmat[x,1]:libmat[x,2],], 2, scale_2)
        rbind(y, NA)
    })
    #tsscale <- do.call(rbind, tmp)
    tsscale <- apply(ts, 2, scale_2)

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## -- NA index
    index <- 1:nrow(ts); names(index) <- rownames(ts)
    naindex <- index[rownames(smlsub[smlsub$missing=='N', ])]

    ## ============================================= ##
    smapres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=13
        
        asv <- tsscale[,k]
        absasv <- ts[,k]

            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## -- Optimaization parameter
            E <- bestE_bidirect(asv, E_range = E.range, lib = libmat, criteria = "rmse", show_fig = F)
            
            ts_tmp <- embed(c(rep(NA,  E-1), asv),  E)
            theta <- bestT_bidirect_block_lnlp(ts_tmp, theta_range=theta_test, lib = libmat, criteria = "rmse", show_fig = F)
               
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
                        
            smapout <- s_map(asv, lib = libmat, pred=libmat, 
                             E = E, theta= theta, tp = tp, tau = 1,
                             stats_only=FALSE)
            
            
            
            stat <- cbind(asv=colnames(tsscale)[k], smapout[,c(3, 1, 5)])
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
        
        
        stat
    }    
    
    ## ============================== ##

    statall <- rbind(statall, cbind(treat=i, smapres))
    
    ## ============================== ##
}


saveRDS(statall, 'Table/05_05_smap_statics')    

## ============================== ##
stattmp <- readRDS('Table/05_05_smap_statics')

thetamat <- as.matrix( table(stattmp[,c(1,5)]) )
a <- 1-thetamat[,1]/rowSums(thetamat)
mean(a)

bardf <- as.data.frame( table(stattmp[,c(1,5)]) )

bardf$theta= factor(bardf$theta, levels= rev(theta_test) )

bardf$treat <- gsub('medium ', 'Medium-', bardf$treat)
bardf$treat <- factor(bardf$treat, levels=rev(c("Water/Medium-A",  "Water/Medium-B","Water/Medium-C",
                                            "Soil/Medium-A", "Soil/Medium-B",   "Soil/Medium-C" )))

g <- ggplot(bardf)+
     geom_bar(aes(x=treat, y= Freq, fill= theta), size=1.4, stat='identity', position='fill')+
     theme_classic(base_size=40)+
     guides(fill=guide_legend(nrow=1,label.position = 'bottom'))+
     theme(#panel.border=element_rect(fill=NA, size=1),
           text=element_text(family='Arial'),
           legend.position='bottom',
           legend.key.width = unit(0.5,'cm'),
           legend.text = element_text(size=10,angle=30, hjust=0.5, vjust=1),
           plot.margin = unit(c(1, 3, 1, 1), "lines"))+
     scale_fill_brewer(palette='Spectral')+
     labs(x='', y='Proportion')+
     #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
     scale_y_continuous(expand=c(0,0))+
     coord_flip()

ggsave(plot=g, filename=sprintf('%s/theta_bar.tiff', dir$figdir),
       w=14, h=7)


## ============================== ##


color <- readRDS('Table/03_color_palette.rds')

g <- ggplot(statall)+
     geom_violin(aes(x=treat, y=theta, color=treat), size=1.4)+
     theme_minimal(base_size=20)+
     theme(panel.border=element_rect(fill=NA, size=1),
           text=element_text(family='Arial'))+
     scale_color_manual(values=color$treat)+
     labs(x='', y='THeta')+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggsave(plot=g, filename=sprintf('%s/theta_distribution.tiff', dir$figdir),
       w=12, h=7)

## ============================== ##
head(statall)

color <- readRDS('Table/03_color_palette.rds')

g <- ggplot(statall)+
    geom_violin(aes(x=treat, y=theta, color=treat), size=1.4)+
    theme_minimal(base_size=20)+
    theme(panel.border=element_rect(fill=NA, size=1),
          text=element_text(family='Arial'))+
    scale_color_manual(values=color$treat)+
    labs(x='', y='THeta')+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggsave(plot=g, filename=sprintf('%s/theta_distribution.tiff', dir$figdir),
       w=7, h=12)


g <- ggplot(statall[statall$tp %in% c(1,7),],aes(x=as.factor(theta), y=rho))+
    geom_violin(size=1.4)+
    stat_summary(size=3, geom='point', fun='mean')+
    theme_minimal(base_size=26)+
    theme(panel.border=element_rect(fill=NA, size=1),
          text=element_text(family='Arial'))+
    scale_color_manual(values=color$treat)+
    labs(x='', y='THeta')+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    facet_wrap(~tp, ncol=1)+
    coord_flip()

ggsave(plot=g, filename=sprintf('%s/theta_error.tiff', dir$figdir),
       w=12, h=18)
