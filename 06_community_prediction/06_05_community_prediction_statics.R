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

load.lib( c('rEDM', 'ggplot2', 'tidyr', 'cowplot', 'scales', 'extrafont', 'RColorBrewer', 'gridExtra', 'vegetarian'))

# -- Create directory to save
dir <- make.dir('06_community_prediction/06_05_output')

# -- Load data table
dlist  <- readRDS("Table/03_matrixList.rds")
libmat=dlist[[7]]
sml <- readRDS('Table/06_04_sample_info.rds')

pred.comm <- readRDS("Table/06_04_predict_community.rds")

color <- readRDS('Table/03_color_palette.rds')
############################################################################
pred.fill=c(simplex='royalblue2', smapT='steelblue1', smap0='firebrick1', AR='orange1', null='darkolivegreen3')
pred.color=c(simplex='royalblue4', smapT='skyblue4', smap0='firebrick3', AR='darkorange4', null='darkolivegreen4')

target <- c('simplex', 'smapT', 'smap0', 'null')
abrupt <- 'abruptness_rel_tw5_tp1'

erts <- erdist <- erdist2 <- erfocus <- erfocus2 <- delta_focus <- summs <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ======================================= ##
    ## -- Extract treatment
    obs <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ]
    
    pred <- c() 
    for(l in 1:10){
        cols <- unlist(sapply(target, function(x) grep(x, colnames(pred.comm[[i]][[1]])) ))
        pred <- rbind(pred, cbind(tp=l, pred.comm[[i]][[l]][,cols] ))
    }
    ## ======================================= ##
    ## -- Error type
    for(n in c('rel')){ #n='rel'
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Extract target column
        errorscol <- c()
        for(l in target){
            
            errorscol <- c(errorscol, grep(sprintf('^%s_%s',l, n), colnames(smlsub)))
        }
        
        errorsml <- smlsub[,c(1:160, errorscol, grep(sprintf('delta_.*_%s', n), colnames(smlsub)) )]
        
        ## -- Abrupt point index
        sml_re <- c()
        for(r in 1:8){ #r=2
            sub <- smlsub[smlsub$replicate.id==r, ]
            sub$abrupt_index <- 1:110-order(sub[,abrupt], decreasing=TRUE)[1]
            sml_re <- rbind(sml_re, sub)
        }
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Visualising error time-series
        lferor <- gather(sml_re, key, value, errorscol)
        lferor2 <- cbind(lferor, group=do.call(rbind, strsplit(as.character(lferor$key), sprintf('_%s_tp_', n))))
        lferor2$group.2 <- factor(lferor2$group.2, levels=c(1:10))
        error_ts <- ggplot(lferor2[lferor2 $group.2%in%c(1,4,7,10),], aes(x=index))+
                    geom_point(aes(y= abruptness_abs_tw5_tp1), size=0.8, color='grey40')+
                    geom_line(aes(y= abruptness_abs_tw5_tp1),color='grey40')+
                    geom_point(aes(y= value, color=group.1), size=0.8)+
                    geom_line(aes( y= value, color= group.1, group= group.1), size=0.5)+
                    facet_wrap(~group.2, scales='free', ncol=1)+
                    scale_color_manual(values=pred.fill[target])+
                    theme_minimal()+
                    theme(panel.border=element_rect(size=0.8, fill=NA),
                          text=element_text(family='Arial'))+
                    labs(title=i, y='Error (euclid distance)', x='Time', caption='Black line represent Abruptness')
        
        ## -- Visualising error distribution
        lferor2$group.1 <- factor(lferor2$group.1, levels=target)
        errorbox <- ggplot(lferor2[which(lferor2$missing=='N'), ], aes(x= group.1, y= value, color= group.1, fill= group.1))+
                    geom_jitter(alpha=0.5)+
                    geom_violin(trim=TRUE, alpha=0.8)+
                    geom_point( stat = "summary", fun = "mean", color='black')+
                    scale_color_manual(values=pred.color[target])+
                    scale_fill_manual(values=pred.fill[target])+
                    facet_wrap(~group.2, nrow=2, dir='h')+
                    theme_minimal()+
                    theme(panel.border=element_rect(size=0.8, fill=NA),
                          text=element_text(family='Arial'))+
                    labs(title=i)
        errorbox2 <- ggplot(lferor2[which(lferor2$missing=='N'), ], aes(x= group.1, y= value))+
                    geom_jitter(alpha=0.5, aes(fill=abrupt_index), shape=21)+
                    #geom_violin(aes(color= group.1), trim=TRUE, alpha=0.8)+
                    geom_point( stat = "summary", fun = "mean", color='black')+
                    scale_color_manual(values=pred.color[target])+
                    scale_fill_gradient2()+
                    facet_wrap(~group.2, nrow=2, dir='h')+
                    theme_minimal()+
                    theme(panel.border=element_rect(size=0.8, fill=NA),
                          text=element_text(family='Arial'))+
                    labs(title=i)
        
        ## -- Community assembly
        tmp <- cbind(sml_re, obs[rownames(sml_re), ])
        lfass <- gather(tmp, key, value, -c(1:ncol(sml_re)))
        lfass$key <- factor(lfass$key, level=colnames(obs)[ order(colSums(obs[rownames(sml_re), ]), decreasing=TRUE) ])
        
        assembly <- ggplot(lfass)+
                    geom_bar(aes(x=as.numeric(time), y=value, fill=key), color='grey30', size=0.2,
                             position='fill', stat='identity', show.legend = FALSE, width=1)+
                    geom_point(data=lfass[which(lfass$abrupt_index==0),], aes(x=as.numeric(time), y=1.1), fill='grey90', shape=25)+
                    facet_wrap(~replicate.id, ncol=1)+
                    scale_fill_manual(values=color$asv$ID)+
                    theme_minimal()+
                    theme(text=element_text(family='Arial'))+
                    labs(title=i, y='', x='Time', caption='triangle represent maximum abruptness')
        
        ## -- Focus on Maximum abrupt point
        index_ts <- ggplot(lferor2[lferor2 $group.2%in%c(1,4,7,10),], aes(x=abrupt_index))+
                    geom_point(aes(y= abruptness_abs_tw5_tp1), size=0.8, color='grey40')+
                    geom_line(aes(y= abruptness_abs_tw5_tp1),color='grey40')+
                    geom_point(aes(y= value, color=group.1), size=0.8)+
                    geom_line(aes( y= value, color= group.1, group= group.1), size=0.5)+
                    facet_grid(replicate.id~group.2, scales='free')+
                    scale_color_manual(values=pred.fill[target])+
                    theme_minimal()+
                    theme(panel.border=element_rect(size=0.8, fill=NA),
                          text=element_text(family='Arial'))+
                    labs(title=i, y='Error (euclid distance)', x='Time', color='method',caption='Black line represent Abruptness')+
                    xlim(-10,10)
        
        ## -- Focus on Maximum abrupt point 2
        BoxMeanQuant <- function(x) {
            v <- c(min(x), quantile(x, 0.25), median(x), quantile(x, 0.75), max(x))
            names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
            v
        }
        index_ts2 <- ggplot(lferor2[lferor2$group.2%in%c(1,7) & lferor2$abrupt_index%in%c(-5:5),],
                           aes(x=as.factor(abrupt_index), y= value, color=group.1))+
                     stat_summary(fun.data = BoxMeanQuant, geom = "boxplot",width=0.6, position=position_dodge(width = .7), size=0.6)+
                     #geom_boxplot(aes(x=as.factor(abrupt_index), y= value, color=group.1), size=0.8)+
                     facet_wrap(~group.2, scales='free', ncol=1)+
                     scale_color_manual(values=pred.fill[target])+
                     theme_bw(base_size=15)+
                     theme(panel.border=element_rect(size=0.8, fill=NA),
                           text=element_text(family='Arial'),
                           strip.background = element_blank(),
                           legend.position='bottom')+
                     labs(subtitle=i, y='Error (Bray-Curtis dissimilarity)', x='Time step from the maximal change', color='method')
        
        erts[[sprintf('%s_%s', i, n)]] <-error_ts
        erdist[[sprintf('%s_%s', i, n)]] <- errorbox
        erdist2[[sprintf('%s_%s', i, n)]] <- errorbox2
        erfocus[[sprintf('%s_%s', i, n)]] <- plot_grid(assembly, index_ts, ncol=2, rel_widths=c(0.6, 1))
        erfocus2[[sprintf('%s_%s', i, n)]] <- index_ts2
        
        summs <- rbind(summs, sml_re)
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Visualising difference EDM vs. null

        lfdelta <- gather(sml_re, key, value, grep(sprintf('delta_.*_%s', n), colnames(sml_re)) )
        
        lfdelta2 <- cbind(lfdelta, group=do.call(rbind, strsplit(as.character(lfdelta$key), sprintf('_error_%s_tp_', n))) )
        lfdelta2$group.2 <- factor(lfdelta2$group.2, levels=c(1:10))
        lfdelta2$group.1 <- gsub('delta_', '', lfdelta2$group.1)
        lfdelta2$group.1 <- gsub('smap', 'smapT', lfdelta2$group.1)
        
        ## -- Focus on Maximum abrupt point
        index_ts <- ggplot(lfdelta2[lfdelta2 $group.2%in%c(1,4,7,10),], aes(x=abrupt_index))+
                    geom_hline(yintercept=0)+geom_vline(xintercept=0, linetype=2, alpha=0.7)+
                    geom_point(aes(y= abruptness_abs_tw5_tp1), size=0.8, color='grey40')+
                    geom_line(aes(y= abruptness_abs_tw5_tp1),color='grey40')+
                    geom_point(aes(y= value, color=group.1), size=0.8)+
                    geom_line(aes( y= value, color= group.1, group= group.1), size=0.5)+
                    facet_grid(replicate.id~group.2, scales='free')+
                    scale_color_manual(values=pred.fill[target])+
                    theme_minimal()+
                    theme(panel.border=element_rect(size=0.8, fill=NA),
                          text=element_text(family='Arial'))+
                    labs(title=i, y='Error (Bray-Curtis dissimilarity)', x='Time step from the maximal change', color='method',caption='Black line represent Abruptness')+
                    xlim(-10,10)+
                    labs(title=i)
        
        delta_focus[[sprintf('%s_%s', i, n)]] <- index_ts
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

    }
    
}

## ===================================================== ##

pdf(sprintf('%s/error_timeSeries.pdf', dir$figdir), w=17)
for(l in 1:length(erts)){ plot(erts[[l]])}; dev.off()    

pdf(sprintf('%s/error_distribution.pdf', dir$figdir), w=17)
for(l in 1:length(erdist)){ plot(erdist[[l]])}; dev.off()   

pdf(sprintf('%s/error_distribution_ver2.pdf', dir$figdir), w=17)
for(l in 1:length(erdist)){ plot(erdist2[[l]])}; dev.off()   

pdf(sprintf('%s/error_around_abrupt.pdf', dir$figdir), h=13, w=15)
for(l in 1:length(erfocus)){ plot(erfocus[[l]])}; dev.off()   

pdf(sprintf('%s/error_around_abrupt_sum.pdf', dir$figdir), h=6, w=4)
for(l in 1:length(erfocus2)){ plot(erfocus2[[l]])}; dev.off()  

tiff(sprintf('%s/error_around_abrupt_sum.tiff', dir$figdir), h=7, w=5)
for(l in 1:length(erfocus2)){ plot(erfocus2[[l]])}; dev.off()  

pdf(sprintf('%s/delta_error_timeSeries.pdf', dir$figdir), w=8)
for(l in 1:length(delta_focus)){ plot(delta_focus[[l]])}; dev.off()   


## ===================================================== ##

lfdelta <- gather(summs[summs$missing=='N', ], key, value, grep(sprintf('delta_.*_%s', n), colnames(sml_re)) )

lfdelta2 <- cbind(lfdelta, group=do.call(rbind, strsplit(as.character(lfdelta$key), sprintf('_error_%s_tp_', n))) )
lfdelta2$group.2 <- factor(lfdelta2$group.2, levels=c(1:10))
lfdelta2$group.1 <- gsub('delta_', '', lfdelta2$group.1)
lfdelta2$group.1 <- gsub('smap', 'smapT', lfdelta2$group.1)

pdf(sprintf('%s/nmds_abs.pdf', dir$figdir), w=11)
for(t in 1:10){ #t=1
    delta_mdsabs <- ggplot(lfdelta2[lfdelta2$group.2==t & as.numeric(lfdelta2$time)%in%c(21:110), ])+
        geom_point(aes(x=abs.MDS1, y=abs.MDS2, color=value, shape=as.factor(replicate.id)), size=1.8 )+
        #geom_point(aes(x=abs.MDS1, y=abs.MDS2, fill=value, shape=as.factor(replicate.id)))+
        facet_wrap(~treat1, scales='free')+
        scale_shape_manual(values=c(16:18, 15, 3, 7:9))+
        scale_color_gradientn(colors=c('darkorange3', 'white','darkblue'), na.value=NA, 
                             values=rescale(c(min(lfdelta2$value, na.rm=TRUE), 0, max(lfdelta2$value, na.rm=TRUE))))+
        scale_fill_gradientn(colors=c('darkorange3', 'white','darkblue'), na.value=NA, 
                             values=rescale(c(min(lfdelta2$value, na.rm=TRUE), 0, max(lfdelta2$value, na.rm=TRUE))))+
        theme_minimal()+
        theme(panel.border = element_rect(fill=NA, size=0.8))+
        labs(shape='ID', subtitle=sprintf('tp=%s',t))
        
        
    plot(delta_mdsabs)
}
dev.off()

pdf(sprintf('%s/nmds_rel.pdf', dir$figdir), w=11)
for(t in 1:10){ #t=1
    delta_mdsabs <- ggplot(lfdelta2[lfdelta2$group.2==t & as.numeric(lfdelta2$time)%in%c(21:110), ])+
        geom_point(aes(x=rel.MDS1, y=rel.MDS2, color=value, shape=as.factor(replicate.id)), size=1.8 )+
        #geom_point(aes(x=abs.MDS1, y=abs.MDS2, fill=value, shape=as.factor(replicate.id)))+
        facet_wrap(~treat1, scales='free')+
        scale_shape_manual(values=c(16:18, 15, 3, 7:9))+
        scale_color_gradientn(colors=c('darkorange3', 'white','darkblue'), na.value=NA, 
                              values=rescale(c(min(lfdelta2$value, na.rm=TRUE), 0, max(lfdelta2$value, na.rm=TRUE))))+
        scale_fill_gradientn(colors=c('darkorange3', 'white','darkblue'), na.value=NA, 
                             values=rescale(c(min(lfdelta2$value, na.rm=TRUE), 0, max(lfdelta2$value, na.rm=TRUE))))+
        theme_minimal()+
        theme(panel.border = element_rect(fill=NA, size=0.8))+
        labs(shape='ID', subtitle=sprintf('tp=%s',t))
    
    
    plot(delta_mdsabs)
}
dev.off()

## ===================================================== ##

ela <- readRDS( 'Table/04_03_sample_info.rds')
head(ela)

df <- cbind(sml, ela[rownames(sml), -13:0+ncol(ela)])

head(sml)
g <- ggplot(df)+
	 geom_point(aes(x= bray_against_replicate, y=smapT_abs_tp_1, color=SS.Entropy))+
	 scale_color_gradientn( colors=brewer.pal(11,'Spectral'))+
	 facet_wrap(~treat1)+
	 theme_minimal(base_size=25)+
	 theme(panel.border = element_rect(fill=NA, size=0.8),
	 	   text=element_text(family='Arial'))+
	 labs(x='Beta diversity against replicate', y='Prediction error (S-map)')
	 
ggsave(plot=g, filename=sprintf('%s/betaDiversity_error.tiff', dir$figdir),
	   w=12, h=7)
	   
	   
## ===================================================== ##

lf <- gather(df, key, value, which(colnames(df) %in% c('smapT_abs_tp_1', 'smapT_abs_tp_7') ))
	   
g <- ggplot(lf, aes(x=as.factor(gsub('smapT_abs_tp_', '', key)), y=value, color=treat1))+
	 stat_summary(fun.data = BoxMeanQuant, geom = "boxplot",width=0.6, position=position_dodge(width = .7), size=0.9)+
	 scale_color_manual(values=color$treat)+
	 theme_minimal(base_size=25)+
	 theme(panel.border = element_rect(fill=NA, size=0.8),
	 	   text=element_text(family='Arial'))+
	 labs(x='Time steps in forecasting', y='Prediction error (S-map)')
ggsave(plot=g, filename=sprintf('%s/Error_distribution.tiff', dir$figdir),
	   w=12, h=7)
	   
g <- ggplot(lf, aes(x=as.factor(gsub('smapT_abs_tp_', '', key)), y=value, color=treat1))+
	 geom_violin(size=1.2)+
	 stat_summary(geom='point', fun='mean', position=position_dodge(width = .9), size=3)+
	 scale_color_manual(values=color$treat)+
	 theme_minimal(base_size=25)+
	 theme(panel.border = element_rect(fill=NA, size=0.8),
	 	   text=element_text(family='Arial'))+
	 labs(x='Time steps in forecasting', y='Prediction error (S-map)')
ggsave(plot=g, filename=sprintf('%s/Error_distribution_ver.Violin_.tiff', dir$figdir),
	   w=12, h=4)	   	  	   


## ===================================================== ##
## GLM and ANOVA 
library(lme4)

stat.all <- cormat <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[2]
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
    
    ts <- na.omit(dlist[[i]]) 
    smlsub <- sml[rownames(ts),]
    smlsub <- smlsub[order(smlsub$time), ]
    
    smlsub$shanon <- apply(ts, 1, d, lev = "alpha", q = 1)
    
    stat.all <- rbind(stat.all, smlsub)
}   

library(glmmML)
library(nlme)

stat.all.tmp <- stat.all
stat.all.tmp$replicate.id <- as.factor(stat.all.tmp$replicate.id)

sink(sprintf('%s/community_prediction_anova_main.txt', dir$tabledir))
a <- lmerTest::lmer(smapT_rel_tp_1 ~ as.factor(inoculum) + as.factor(resource) + shanon + bray_against_replicate + (1| replicate.id), data= stat.all.tmp[stat.all.tmp $missing=='N',])
anova(a)
sink()
