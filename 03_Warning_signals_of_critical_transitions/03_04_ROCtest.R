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

load.lib( c('MASS', 'lme4', 'ggplot2', 'tidyr', 'cowplot', 'scales', 'car', 'RColorBrewer', 'extrafont','svglite','ggstar'))

# -- Create directory to save
dir <- make.dir('03_Warning_signals_of_critical_transitions/ROC_Test')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/sample_info.rds')
div <- readRDS('Table/01_Diversity.rds')
ela <- readRDS('Table/ELA.rds')
stability <- readRDS('Table/03_Stability.rds')
#ews <- readRDS('Table/03_Early_warning_signals.rds')
abruptness <- readRDS('Table/01_Abruptness.rds')

color <- readRDS('Table/color_palette.rds')

# -- parameter
abrupttw = 'abruptness_rel_tw5_tp1'
delay = 11

for.parallel(8)

############################################################################
## -- Possible candidate for early warnings signals
ews <- c('Energy.SS.Energy', 'SS.Entropy', 'DS_0710_tp1', 'LSS_0710_tp1')

smlselect <- cbind(abruptness[rownames(sml),'abruptness_rel_tw5_tp1'], 
                   ela[rownames(sml), colnames(ela)%in%ews],
                   stability[rownames(sml), colnames(stability)%in%ews])
colnames(smlselect) <- c('Community change',
                 "Energy gap", "Stable-state entropy",
                 'Local Lyapunov stability','Local structural stability')
                   
smlews <- cbind(sml, smlselect[rownames(sml), ])
smlews <- smlews[which(smlews$missing=='N'),]

plot(sub$`Stable-state entropy`, sub$`Community change`)
## ======================================================= ##

glistal <-c()
for(i in names(dlist)[-c(5:7)]){ # i=names(dlist)[2]
    
    ## =================================================== ##
    ## -- Extract one treatment 
    sub <- smlews[which(smlews$treat1==i), ]

    ## =================================================== ##
    ## -- Time delayed abruptness
    dly <- c()
    for(r in 1:8){
        	
        	smlrep <- sub[which(sub $replicate.id==r),]
        	smlemb <- embed(c(smlrep[,'Community change'], rep(NA,7)), 8)
        	dimnames(smlemb) <- list(rownames(smlrep), paste('abruptness_delay', 7:0, sep='') )
        	dly <- rbind(dly, smlemb)
    }
   
    ## =================================================== ##    
    glist <- c()
    
    ## -- Correlation between abruptness at tp7 and EWS 
    for(d in 1 ){ #d=1
        
        tmp <- sub
        tmp[,'Community change'] <- dly[rownames(tmp),d]
        lf <- gather(tmp, key, value, -c(1:9))
        lf$key <- factor(lf$key, levels=c(
                 "Energy gap", "Stable-state entropy",
                 'Local Lyapunov stability','Local structural stability'
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
            geom_hline(yintercept=0.5, linetype=4, size=0.3)+
            geom_point(aes(y=`Community change`, x=value, color=as.factor(replicate.id)), 
            			show.legend=FALSE, size=0.05)+
            geom_smooth(aes(y= `Community change`, x=value, group=as.factor(replicate.id), 
            			color=as.factor(replicate.id), linetype=as.factor(sig)),
                        se=FALSE, method='lm', show.legend=FALSE, size=0.35)+
            facet_wrap(~key, scales='free', ncol=1,
                       strip.position = "bottom")+
            labs(subtitle=i )+
            scale_color_manual(values=color$replicate)+
            theme_minimal(base_size=6)+
            labs(x=NULL, y='Observed community change (abruptness score)')+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid=element_blank(),
                  strip.text = element_text(size=5, vjust=1, margin=margin(-0.001,0.01,0,0, unit='inch') ),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=5),
                  axis.text.y = element_text(size=5, hjust=1, margin=margin(0,0,-0.05,0, unit='inch')),
	                panel.border=element_rect(fill=NA, size=0.4),
	                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt")  )+
             ylim(c(0,1))
        glist[[sprintf('%s', d)]] <- g
        glistal[[i]] <- g
        
        
        if(i ==names(dlist)[2]){
        	g <- ggplot(lf)+
            geom_hline(yintercept=0.5, linetype=4, size=0.3)+
            geom_point(aes(y=`Community change`, x=value, color=as.factor(replicate.id)), show.legend=FALSE, size=0.05)+
            geom_smooth(aes(y= `Community change`, x=value, group=as.factor(replicate.id), 
            			color=as.factor(replicate.id), linetype=as.factor(sig)),
                        se=FALSE, method='lm', show.legend=FALSE, size=0.35)+
            facet_wrap(~key, scales='free', nrow=2,
                       strip.position = "bottom")+
            labs(subtitle=i )+
            scale_color_manual(values=color$replicate)+
            theme_minimal(base_size=8)+
            labs(x=NULL, y='Observed community change (abruptness score)', subtitle='7-day-ahead forecasting results (Soil/Medium-A)')+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid=element_blank(),
                  strip.text = element_text(size=7, vjust=1, margin=margin(-0.001,0.01,0,0, unit='inch') ),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=6),
                  axis.text.y = element_text(size=6, hjust=1, margin=margin(0,0,-0.05,0, unit='inch')) ,
	                panel.border=element_rect(fill=NA, size=0.4),
	                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt") )+
        	scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
             ylim(c(0,1))
        	ggsave(filename = sprintf('%s/abruptness_indices_SA.tiff', dir$figdir, gsub('/', '_', i)),
       		plot=g, h=9, w=12)

			svglite(file=sprintf('%s/abruptness_indices_SA.svgz', dir$figdir),
	   		 height =8/2.5, width = 8/2.5)
			plot(g)
			dev.off()
        }  
		
    }

}

ggsave(filename = sprintf('%s/abruptness_indices_extend.tiff', dir$figdir, gsub('/', '_', i)),
       plot=plot_grid(plotlist=glistal, nrow=1), h=10, w=10)
svglite(file=sprintf('%s/abruptness_indices_merge.svgz', dir$figdir),
	   		 height =12/2.5, width = 18/2.5)
			plot(plot_grid(plotlist= glistal, nrow=1))
			dev.off()

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
lists <- lists2  <- c()
for(l in colnames(smlews)[10:ncol(smlews)]){ #l= colnames(smlews)[10:ncol(smlews)][2]
    
    ## =============================================== ##
    ggdf <- ggauc <- exact <- c()
    for( i in names(dlist)[-c(5:7)]){ #i=names(dlist[2])
        
        smlsub <- smlews[which(smlews$treat1==i),]
        
        ## =========================================== ##
   		## -- Time delayed abruptness
        dly <- c()
        for(r in 1:8){
        	
        	smlrep <- smlsub[which(smlsub$replicate.id==r),]
        	smlemb <- embed(c(smlrep[,'Community change'], rep(NA,7)), 8)
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
    ggdf$treat <- factor(ggdf$treat, levels=c('Water/Medium-A', 'Soil/Medium-A',  'Water/Medium-B', 'Soil/Medium-B'))
    ggauc$treat <- factor(ggauc$treat, levels=c('Water/Medium-A', 'Soil/Medium-A',  'Water/Medium-B', 'Soil/Medium-B'))
    g <- ggplot(ggdf, aes(x=1-specify, y=sense))+
        geom_abline(linetype=4, color='grey70', size=0.3)+
        geom_path( size=0.35)+
        geom_star(data=ggdf, aes(y=youden), fill='darkorange',size=2, starshape=1)+
        geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC = %.3f',auc)), size=2 )+
        facet_wrap(~treat, nrow=1)+
        theme_minimal(base_size=8)+
        labs(y='Detection rate', x='False detection rate')+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid=element_blank(),
              strip.text = element_text(size=7, hjust=0),
              strip.text.x = element_text( margin = margin( r=1, l=1, b = 0, t = 0) ) ,
              text=element_text(family='Arial'),
              #axis.text=element_text(size=15),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=6),
              panel.border=element_rect(fill=NA, size=0.6),
                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt"))+
        	scale_x_continuous(guide = guide_axis(check.overlap = TRUE), breaks=seq(0, 1, 0.2))+
        	coord_fixed()
    g3 <- ggplot(ggdf)+
    		geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=0.3, linetype=2)+
            geom_path( size=0.4, aes(y=sense, x=exp), color='blue')+
            geom_path( size=0.4, aes(y=1-specify, x=exp), color='firebrick3')+
            #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
            #geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Threshold=%s', round(youden2, 3))), 
             #         size=8,hjust=0)+
            facet_wrap(~facet, nrow=1)+
            theme_minimal(base_size=8)+
            labs(y='Detection rate or\nFalse detection rate', x=l)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid=element_blank(),
                  strip.text = element_text(hjust=0, size=7),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=6),
                  axis.text.y = element_text(size=6),
              panel.border=element_rect(fill=NA, size=0.4),
                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt"))+
        	scale_x_continuous(guide = guide_axis(check.overlap = TRUE))+
        	coord_fixed()
        	
    if(any(c('DS_0710_tp1', 'LSS_0710_tp1')==l)){
    g3 <- ggplot(ggdf)+
    		geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=0.3, linetype=2)+
            geom_path( size=0.4, aes(y=sense, x=exp), color='blue')+
            geom_path( size=0.4, aes(y=1-specify, x=exp), color='firebrick3')+
            #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
            #geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Optimized\nthreshold=%s', round(youden2, 3))), 
            #          size=8,hjust=0)+
            facet_wrap(~facet, nrow=1, scales='free_x')+
            theme_bw(base_size=8)+
            labs(y='Detection rate/\nFalse detection rate', x=l)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid=element_blank(),
                  strip.text = element_text(size=7,hjust=0),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=6),
                  axis.text.y = element_text(size=6))+
        	scale_x_continuous(guide = guide_axis(check.overlap = TRUE))}        	    	
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	title <- ggdraw() + 
			  draw_label(gsub('\n', ' ', l),
			    x = 0,
			    hjust = 0,
			    size=7
			  ) +
			  theme(
			    # add margin on the left of the drawing canvas,
			    # so title is aligned with left edge of first plot
			    plot.margin = margin(0, 0, 0, 29)
			  )
			  
	title2 <- ggdraw() + 
			  draw_label(gsub('\n', ' ', l),
			    x = 0,
			    hjust = 0,
			    size=7
			  ) +
			  theme(
			    # add margin on the left of the drawing canvas,
			    # so title is aligned with left edge of first plot
			    plot.margin = margin(0, 0, 0, 23)
			  )		  
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##      
	svglite(file=sprintf('%s/roc_AUC_%s.svgz', dir$figdir, gsub('\n', ' ', l)),
	   		 height =3.5/2.5, width = 10/2.5)
	plot(plot_grid(title2, g, ncol=1, rel_heights=c(0.07, 1)))
	dev.off()
	
	svglite(file=sprintf('%s/roc_threshold_%s.svgz', dir$figdir, gsub('\n', ' ', l)),
	   		 height =3.5/2.5, width = 10/2.5)
	plot(plot_grid(title, plot_grid(g3, ncol=1, align='v'), ncol=1, rel_heights=c(0.07, 1)) )
	dev.off()
	
	lists[[l]] <-plot_grid(title2, g, ncol=1, rel_heights=c(0.07, 1))
	lists2[[l]] <-plot_grid(title, g, ncol=1, rel_heights=c(0.07, 1))

}

svglite(file=sprintf('%s/roc_AUC_all.svgz', dir$figdir, gsub('\n', ' ', l)),
	   		 height =15/2.5, width = 10/2.5)
plot(plot_grid(plotlist=lists, ncol=1))
dev.off()
## =================================================== ##
## -- ROC test for Local Lyapunov stability and Stable state entropy in all treatment

for(l in colnames(smlews)[c(12,11)]){ #l=colnames(smlews)[c(12)]
    
    ggdf <- ggauc <- exact <- c()

    smltmp <- smlews
    smlsub <- c()
    for(i in names(dlist)[-c(5:7)]){
    	
    	tmp <- smltmp[which(smltmp $treat1==i),]
        dly <- c()
        for(r in 1:8){
        	
        	smlrep <- tmp[which(tmp $replicate.id==r),]
        	smlemb <- embed(c(smlrep[,'Community change'], rep(NA,7)), 8)
        	dimnames(smlemb) <- list(rownames(smlrep), paste('abruptness_delay', 7:0, sep='') )
        	dly <- rbind(dly, smlemb)
        }
        
    	smlsub <- rbind(smlsub, cbind(tmp, dly[rownames(tmp), ]))
    }
     
    roctmp <- roctest(data=smlsub, exp=l, responce = 'abruptness_delay7', responceth = 0.5)
        
    exact <- rbind(exact, roctmp[[2]])
    ggdf <- rbind(ggdf, cbind(treat=i, roctmp[[1]]))
    ggauc <- rbind(ggauc, data.frame(treat=i,x=0.7, y=0.2, auc=roctmp[[3]]))
    subtitle =sprintf('Threshold = %.3f', na.omit(ggdf)$youden2)
    
    write.csv(exact, sprintf('%s/exact_test_%s.csv', dir$tabledir,l), row.names=FALSE)
    ggdf$treat <- factor(ggdf$treat, levels=c('Water/Medium-A', 'Soil/Medium-A',  'Water/Medium-B', 'Soil/Medium-B'))
    ggauc$treat <- factor(ggauc$treat, levels=c('Water/Medium-A', 'Soil/Medium-A',  'Water/Medium-B', 'Soil/Medium-B'))
    
    g <- ggplot(ggdf, aes(x=1-specify, y=sense))+
        geom_abline(linetype=4, color='grey70', size=0.3)+
        geom_path( size=0.35)+
        geom_star(data=ggdf, aes(y=youden), fill='darkorange',size=2, starshape=1)+
        geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC = %.3f',auc)), size=2.3)+
        theme_minimal(base_size=7)+
        labs(y='Detection rate', x='False detection rate')+
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid=element_blank(),
              strip.text = element_text(size=6, vjust=1 ),
              text=element_text(family='Arial'),
              #axis.text=element_text(size=15),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=6) ,
              panel.border=element_rect(fill=NA, size=0.4),
                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt") )+
        	scale_x_continuous(guide = guide_axis(check.overlap = TRUE), breaks=seq(0, 1, 0.2))
      g3 <- ggplot(ggdf)+
    		geom_vline(data=ggdf, aes(xintercept=youden2), color='orangered',size=0.3, linetype=2)+
            geom_path( size=0.4, aes(y=sense, x=exp), color='blue')+
            geom_path( size=0.4, aes(y=1-specify, x=exp), color='firebrick3')+
            geom_star(data=ggdf, aes(x=youden2,y=youden), fill='darkorange',size=2, starshape=1)+
            #geom_text(data=ggauc,aes(x=x, y=y, label=sprintf('AUC=%s',round(auc, 3))), size=8 )+
            geom_text(data=na.omit(ggdf),aes(x=0, y=0.2, label=sprintf('Threshold\n= %s', round(youden2, 3))), 
                     size=2.3,hjust=0)+
            #facet_wrap(~facet, nrow=1)+
            theme_minimal(base_size=7)+
            labs(y='Detection rate or\nFalse detection rate', x=l)+
            theme(strip.background = element_blank(),
                  strip.placement = "outside",
                  panel.grid=element_blank(),
                  strip.text = element_text(hjust=0, size=6),
                  text=element_text(family='Arial'),
                  #axis.text=element_text(size=15),
                  axis.text.x = element_text(size=6),
                  axis.text.y = element_text(size=6) ,
	                panel.border=element_rect(fill=NA, size=0.4),
	                axis.ticks.x = element_line(size = 0.2),
                   		axis.ticks.y = element_line(size = 0.2),
                   	axis.ticks.length = unit(1.8, "pt") )+
        	xlim( c(0, ggdf$exp[which.min(ggdf$sense)]+0.5 ) )
        
     quartz(type="pdf", file=sprintf('%s/roc_test2_%s_all.pdf', dir$figdir, l),
	   		 height = 4/2.5, width = 10/2.5, pointsize=6)
	plot(plot_grid(g,g3, ncol=2, align='hv'))
	dev.off()
	
	svglite(file=sprintf('%s/roc_test2_%s_all.svgz', dir$figdir, l),
	   		 height =3/2.5, width = 8/2.5)
	plot(plot_grid(g,g3, ncol=2, align='hv'))
	dev.off()}

