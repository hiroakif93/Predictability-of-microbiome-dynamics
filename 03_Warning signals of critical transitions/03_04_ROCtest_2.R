############################################################################
####
#### R script for Fujita (2019)
####
#### Finding expected early warnings signal by visualizing
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')

load.lib( c('MASS', 'lme4', 'ggplot2', 'tidyr', 'cowplot', 'scales', 'car', 'RColorBrewer', 'extrafont'))

# -- Create directory to save
dir <- make.dir('03_Warning signals of critical transitions/ROC_Test')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/sample_info.rds')
div <- readRDS('Table/01_Diversity.rds')
ela <- readRDS('Table/ELA.rds')
stability <- readRDS('Table/03_Stability.rds')
ews <- readRDS('Table/03_Early_warning_signals.rds')
abruptness <- readRDS('Table/01_Abruptness.rds')

color <- readRDS('Table/color_palette.rds')

sml <- unique(cbind(sml, div[rownames(sml),],  ews[rownames(sml),], stability[rownames(sml), c('DS_0710_tp1', 'LSS_0710_tp1')], ela[rownames(sml),]))

# -- parameter
abrupttw = 'abruptness_rel_tw5_tp1'
delay = 11

for.parallel(8)

############################################################################
## -- Possible candidate for early warnings signals
ews <- c('cv_tp10', 'skewness_tp10', 'DS_0710_tp1', 'LSS_0710_tp1', 'Energy.SS.Energy', 'SS.Entropy', 'shannon')

smlselect <- cbind(abruptness[rownames(sml),'abruptness_rel_tw5_tp1'], sml[, ews])
colnames(smlselect) <- c('Community change',
                 "Mean coefficient of variation\nof ASV abundance",'Mean skewness of\nASV abundance',
                 'Local Lyapunov stability','Local structural stability',
                 "Energy gap", "Stable state entropy",
                  "Shannon's Entropy")
                   
smlews <- cbind(sml[,1:10], smlselect[rownames(sml), ])

## ======================================================= ##

glistal <-c()
for(i in names(dlist)[-c(5,7)]){ # i=names(dlist)[2]
    
    ## =================================================== ##
    ## -- Extract one treatment 
    ts <- na.omit(dlist[[i]])
    sub <- smlews[rownames(ts), ]

    ## =================================================== ##
    ## -- Time delayed abruptness
    dly <- c()
    for(r in 1:8){
        	
        	smlrep <- sub[which(sub $replicate.id==r),]
        	smlemb <- embed(c(smlrep[,11], rep(NA,7)), 8)
        	dimnames(smlemb) <- list(rownames(smlrep), paste('abruptness_delay', 7:0, sep='') )
        	dly <- rbind(dly, smlemb)
    }
   
    ## =================================================== ##    
    glist <- c()
    
    ## -- Correlation between abruptness at tp7 and EWS 
    for(d in 1 ){ #d=1
        
        tmp <- sub
        tmp[,11] <- dly[rownames(tmp),d]
        lf <- gather(tmp, key, value, -c(1:11))
        lf$key <- factor(lf$key, levels=c('Community change',"Shannon's Entropy",
                 'Mean skewness of\nASV abundance',"Mean coefficient of variation\nof ASV abundance",
                 'Local Lyapunov stability','Local structural stability',
                 "Energy gap", "Stable state entropy"
                  ))
        lf$sig <- NA
        for(v in unique(as.character(lf$key))){ #v=unique(as.character(lf$key))[1]
            
            pvals <- c()
            for(r in 1:8){ #r=1
                lfsub <- lf[ which(lf$key==v & lf$replicate.id==r), ]
                pval <- summary(aov(`Community change`~value, data=lfsub))[[1]]$`Pr(>F)`[1]
                pvals <- c(pvals, pval)
            }
            
            fdr <- p.adjust(pvals, method='fdr')
            for(r in 1:8){ #r=1
           
                lf$sig[which(lf$key==v & lf$replicate.id==r)] <- ifelse(fdr[r] < 0.05, 1, 2)
                
            }
            
        }
        
        g <- ggplot(lf)+
            geom_hline(yintercept=0.5, linetype=4, size=1)+
            geom_point(aes(y=`Community change`, x=value, color=as.factor(replicate.id)), show.legend=FALSE)+
            geom_smooth(aes(y= `Community change`, x=value, group=as.factor(replicate.id), color=as.factor(replicate.id), linetype=as.factor(sig)),
                        se=FALSE, method='lm', show.legend=FALSE)+
            facet_wrap(~key, scales='free', nrow=2,
                       strip.position = "bottom")+
            labs(subtitle=i )+
            scale_color_manual(values=color$replicate)+
            theme_bw(base_size=22)+
            labs(x=NULL, y=NULL)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_line(color='grey90', size=0.3),
                  strip.text = element_text(vjust=1, margin=margin(-0.001,0,0,0.15, unit='inch') ),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=15),
                  axis.text.y = element_text(size=15, hjust=1, margin=margin(0,0,-0.05,0, unit='inch')) )+
             ylim(c(0,1))
        glist[[sprintf('%s', d)]] <- g
        glistal[[i]] <- g
    }

}

ggsave(filename = sprintf('%s/abruptness_indices_extend.tiff', dir$figdir, gsub('/', '_', i)),
       plot=plot_grid(plotlist=glistal, nrow=1), h=8, w=15)

###################################################################################################
## -- ROC test
roctest <- function(data=NULL, exp=1, responce=2, responceth=0){
    
    require(pROC)
    
    ## ----------------------------- ##
    ## -- ROC test
    test <- data[,c(exp, responce)]
    colnames(test) <- c('expvar', 'responcevar')
    test$change <- ifelse( test[,2]>=responceth, 1,0)
    ROC <- roc(change ~ expvar, data = test, ci = TRUE)
    
    ## ----------------------------- ##
    ## -- Younden index
    youden=which.max(ROC$specificities+(ROC$sensitivities-1))
    
    df <- data.frame(specify=ROC$specificities, sense=ROC$sensitivities, 
                     exp=c(NA,sort(unique(na.omit(test)[,'expvar'] ))))
    
    df$youden <- df$youden2 <- NA
    df$youden[youden] <- df[youden,'sense']
    df$youden2[youden] <- df[youden,'exp']
    
    ## ----------------------------- ##
    ## -- Exact test
    mat <- cbind(abrupt=table(test$change),
                 ss.ent=table(ifelse( test[,'expvar'] >=c(NA,sort(unique(na.omit(test)[,'expvar'])))[youden], 1,0)))
    
    excat <- data.frame(treat=i, fisher=fisher.test(mat)$p.value, chi=chisq.test(mat)$p.value)
    ## ----------------------------- ##
    return( list(df,excat, ROC$auc))
}

## =================================================== ##
## -- ROC test in each EWS and each treatment

for(l in colnames(smlews)[12:ncol(smlews)]){ #l= colnames(smlsum)[12:ncol(smlsum)][1]
    
    ## =============================================== ##

    ggdf <- ggauc <- exact <- c()
    for( i in names(dlist)[-c(5:7)]){ #i=names(dlist[2])
        
        smlsub <- smlews[which(smlews $treat1==i),]
        
        ## =========================================== ##
   		## -- Time delayed abruptness
        dly <- c()
        for(r in 1:8){
        	
        	smlrep <- smlsub[which(smlsub$replicate.id==r),]
        	smlemb <- embed(c(smlrep[,11], rep(NA,7)), 8)
        	dimnames(smlemb) <- list(rownames(smlrep), paste('abruptness_delay', 7:0, sep='') )
        	dly <- rbind(dly, smlemb)
        }
        
        smltest <- cbind(smlsub, dly[rownames(smlsub),])
        
        ## =========================================== ##
        roctmp <- roctest(smltest, exp=l, responce = 'abruptness_delay7', responceth = 0.5)        
        
        exact <- rbind(exact, roctmp[[2]])
        ggdf <- rbind(ggdf, cbind(treat=i,facet=sprintf('Threshold = %.3f', na.omit(roctmp[[1]])$youden2), 
                                  roctmp[[1]]))
        ggauc <- rbind(ggauc, data.frame(treat=i,x=0.7, y=0.2, auc=roctmp[[3]]))
    }
    
    write.csv(exact, sprintf('%s/exact_test_%s.csv', dir$tabledir,l), row.names=FALSE)
    
    ## =============================================== ##
	## -- Visualize
    g <- ggplot(ggdf, aes(x=1-specify, y=sense))+
        geom_abline(linetype=4, color='grey70', size=0.8)+
        geom_path( size=1)+
        geom_point(data=ggdf, aes(y=youden), color='firebrick3',size=3)+
        geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC = %.3f',auc)), size=8 )+
        facet_wrap(~treat, nrow=1)+
        theme_bw(base_size=28)+
        labs(y='Detection rate', x='False detection rate')+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color='grey90', size=0.3),
              strip.text = element_text(size=30, vjust=1, margin=margin(-0.001,0,0,0.15, unit='inch') ),
              text=element_text(family='Arial'),
              #axis.text=element_text(size=15),
              axis.text.x = element_text(size=21),
              axis.text.y = element_text(size=21))
    g3 <- ggplot(ggdf)+
    		geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=1.5, linetype=2)+
            geom_path( size=1, aes(y=sense, x=exp), color='blue')+
            geom_path( size=1, aes(y=1-specify, x=exp), color='firebrick3')+
            #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
            #geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Threshold=%s', round(youden2, 3))), 
             #         size=8,hjust=0)+
            facet_wrap(~facet, nrow=1)+
            theme_bw(base_size=28)+
            labs(y='Detection rate/\nFalse detection rate', x=l)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_line(color='grey90', size=0.3),
                  strip.text = element_text(hjust=0, size=24),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=21),
                  axis.text.y = element_text(size=21))
    if(l==colnames(smlews)[c(14)] | l==colnames(smlews)[c(15)]){
    g3 <- ggplot(ggdf)+
    		geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=1.5, linetype=2)+
            geom_path( size=1, aes(y=sense, x=exp), color='blue')+
            geom_path( size=1, aes(y=1-specify, x=exp), color='firebrick3')+
            #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
            #geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Optimized\nthreshold=%s', round(youden2, 3))), 
            #          size=8,hjust=0)+
            facet_wrap(~facet, nrow=1, scales='free_x')+
            theme_bw(base_size=28)+
            labs(y='Detection rate/\nFalse detection rate', x=l)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid.minor=element_blank(),
                  panel.grid.major=element_line(color='grey90', size=0.3),
                  strip.text = element_text(hjust=0),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=21),
                  axis.text.y = element_text(size=21))}    
                  
    ggsave(filename = sprintf('%s/roc_test_%s.pdf', dir$figdir, gsub('\n', ' ', l) ),
           plot=plot_grid(g,g3, ncol=1, align='v'), h=9, w=17 )

}

## =================================================== ##
## -- ROC test for Local Lyapunov stability and Stable state entropy in all treatment

for(l in colnames(smlews)[c(14,17)]){ #l=colnames(smlews)[c(17)]
    
    ggdf <- ggauc <- exact <- c()

    smltmp <- smlews
    smlsub <- c()
    for(i in names(dlist)[-c(5:7)]){
    	
    	tmp <- smltmp[which(smltmp $treat1==i),]
        dly <- c()
        for(r in 1:8){
        	
        	smlrep <- tmp[which(tmp $replicate.id==r),]
        	smlemb <- embed(c(smlrep[,11], rep(NA,7)), 8)
        	dimnames(smlemb) <- list(rownames(smlrep), paste('abruptness_delay', 7:0, sep='') )
        	dly <- rbind(dly, smlemb)
        }
        
    	smlsub <- rbind(smlsub, cbind(tmp, dly[rownames(tmp), ]))
    }
     
    plot(smlsub[,'Stable state entropy'], smlsub[,'abruptness_delay7'])    
    roctmp <- roctest(data=smlsub, exp=l, responce = 'abruptness_delay7', responceth = 0.5)
        
    exact <- rbind(exact, roctmp[[2]])
    ggdf <- rbind(ggdf, cbind(treat=i, roctmp[[1]]))
    ggauc <- rbind(ggauc, data.frame(treat=i,x=0.7, y=0.2, auc=roctmp[[3]]))
    subtitle =sprintf('Threshold = %.3f', na.omit(ggdf)$youden2)
    
    write.csv(exact, sprintf('%s/exact_test_%s.csv', dir$tabledir,l), row.names=FALSE)
    
    
    g <- ggplot(ggdf, aes(x=1-specify, y=sense))+
        geom_abline(linetype=4, color='grey70', size=1.2)+
        geom_path( size=1.5)+
        geom_point(data=ggdf, aes(y=youden), color='firebrick3',size=3)+
        geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=15 )+
        theme_bw(base_size=35)+
        labs(y='Detection rate', x='False detection rate')+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color='grey90', size=0.3),
              strip.text = element_text(size=30, vjust=1, margin=margin(-0.001,0,0,0.15, unit='inch') ),
              text=element_text(family='Arial'),
              #axis.text=element_text(size=15),
              axis.text.x = element_text(size=27),
              axis.text.y = element_text(size=27))
    g3 <- ggplot(ggdf)+
        geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=1.2, linetype=2)+
        geom_path( size=1.5, aes(y=sense, x=exp), color='blue')+
        geom_path( size=1.5, aes(y=1-specify, x=exp), color='firebrick3')+
        #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
        geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Threshold\n= %s', round(youden2, 3))), 
                  size=13,hjust=0)+
        theme_bw(base_size=35)+
        labs(y='Detection rate/\nFalse detection rate', x=l)+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor=element_blank(),
              panel.grid.major=element_line(color='grey90', size=0.3),
              strip.text = element_blank( ),
              text=element_text(family='Arial'),
              #axis.text=element_text(size=15),
              axis.text.x = element_text(size=27),
              axis.text.y = element_text(size=27))+
          xlim(c(0, ggdf$exp[which.min(ggdf$sense)]+0.5))
        
    
    
    ggsave(filename = sprintf('%s/roc_test2_%s_all.pdf', dir$figdir, l),
           plot=plot_grid(g,g3, ncol=1, align='v'), h=14, w=9 )
    ggsave(filename = sprintf('%s/roc_test2_%s_all.tiff', dir$figdir, l),
           plot=plot_grid(g,g3, ncol=1, align='v'), h=14, w=9 )
}
