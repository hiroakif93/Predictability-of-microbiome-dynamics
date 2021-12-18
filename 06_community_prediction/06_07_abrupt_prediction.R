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

load.lib( c('rEDM', 'ggplot2', 'tidyr', 'cowplot', 'scales', 'extrafont', 'RColorBrewer', 'ggpubr', 'ggsignif'))

# -- Create directory to save
dir <- make.dir('06_community_prediction/06_07_output')

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

erts <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[1]
    
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
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        lferor <- gather(sml_re, key, value, errorscol)
        lferor2 <- cbind(lferor, group=do.call(rbind, strsplit(as.character(lferor$key), sprintf('_%s_tp_', n))))
        lferor2$group.2 <- factor(lferor2$group.2, levels=c(1:10))
        smap <- lferor2[grep('smap', lferor2$group.1), ]
		statp <- c()
		for(t in c(1,7)){
			
			smaptp1 <- smap[smap $group.2==t, ]
	        head(smaptp1)
	        
	        dayp <- c()
	        for(day in c(0:10)){ #day=1
			
				abruptafter <- smaptp1[smaptp1$abrupt_index%in% day,]
				
				if( all(table(as.character(na.omit(abruptafter[,c('value','group.1')])[,2]))>1) ) {
				    #pval <- pairwise.t.test(abruptafter$value, abruptafter$group.1, pool.sd=FALSE, p.adj='none')
				    #pval <- pval$p.value
				    
				    group1= abruptafter[abruptafter$group.1=='smapT',]
				    group2= abruptafter[abruptafter$group.1=='smap0',]
				    pval <- t.test(group1$value, group2$value, var.equal=F)
					pval <- data.frame(smapTmean=mean(group1$value, na.rm=TRUE), smapTsd=sd(group1$value, na.rm=TRUE),
					  smap0mean=mean(group2$value, na.rm=TRUE), smap0sd=sd(group2 $value, na.rm=TRUE), 
					  df=pval$parameter, t=pval$statistic,  pval=pval$p.value)

				}else{
				    pval<- data.frame(smapTmean=NA, smapTsd= NA,
					  				  smap0mean= NA, smap0sd= NA ,
					  				  df= NA, t= NA,  pval= NA)
				}
		    	dayp <- rbind(dayp, pval)
		    }
		    
		    statp <- rbind(statp , cbind(group.2=t, day=c(0:10), dayp, fdr=p.adjust(dayp$pval, method='fdr')))
		}		
		statdf <- as.data.frame(statp)
	    #statdf$fdr <- p.adjust(statdf$dayp, method='fdr')

		write.csv(statdf, sprintf('%s/%s.csv', dir$tabledir, gsub('/','_',i)))
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Visualising error time-series
        
        ## -- Focus on Maximum abrupt point
        BoxMeanQuant <- function(x) {
            v <- c(min(x), quantile(x, 0.25), median(x), quantile(x, 0.75), max(x))
            names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
            v
        }
        smap$group3 <- apply(smap[,c('replicate.id', 'group.2', 'group.1')], 1, paste, collapse='_')
        index_ts2 <- ggplot(smap[smap$group.2%in%c(1,7) & smap$abrupt_index%in%c(0:10),],
                           aes(x=as.factor(abrupt_index), y= value, color=group.1))+
                     stat_summary(fun.data = BoxMeanQuant, geom = "boxplot",width=0.6, position=position_dodge(width = .7), size=0.6)+
                     geom_point(position= position_dodge(0.8))+
                     #geom_boxplot(aes(x=as.factor(abrupt_index), y= value, color=group.1), size=0.8)+
                     geom_text(data=statdf, aes(x=as.factor(day),y=0.8, label=ifelse(fdr <0.05,'*','')), color='black')+
                     facet_wrap(~group.2, scales='free', ncol=1)+
                     scale_color_manual(values=pred.fill[target])+
                     theme_bw(base_size=15)+
                     theme(panel.border=element_rect(size=0.8, fill=NA),
                           text=element_text(family='Arial'),
                           strip.background = element_blank(),
                           legend.position='bottom')+
                     labs(subtitle=i, y='Error (Bray-Curtis dissimilarity)', x='Time step from the maximal change', color='method')
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

		erts[[i]] <- index_ts2
    }
    
}

## ===================================================== ##

pdf(sprintf('%s/error_timeSeries_test2.pdf', dir$figdir), w=8)
for(l in 1:length(erts)){ plot(erts[[l]])}; dev.off()    
