############################################################################
####
#### R script for Fujita (2019)
####
#### Prediction by using simplex projection
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder --　setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('ggplot2', 'tidyr', 'cowplot'))

# -- Create directory to save
dir <- make.dir('02_Predicting_microbiome_dynamics/Population_dynamics')

# -- Load library and functions
load.lib(c('ggplot2', 'cowplot', 'RColorBrewer', 'tidyr', 'extrafont','gridExtra'))

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

simplexpc <- readRDS('Table/02_01_simplexModel.rds')    
simplexs <- readRDS('Table/02_01_simplexStatics.rds')    
smapTpc <- readRDS('Table/02_02_smapTModel.rds')    
smap0pc <- readRDS('Table/02_02_smap0Model.rds') 
smaps <- readRDS('Table/02_02_smap_statics.rds')    
meanpc <- readRDS('Table/02_03_meanModel.rds')    
arpc <-readRDS('Table/02_03_ARModel.rds')    
ars <- readRDS('Table/02_03_AR_statics.rds')    

color=readRDS('Table/color_palette.rds')

############################################################################
## -- Compiling forecasting results

statall <- pvalall <- c()

for(i in names(simplexpc)){ #i=names(simplexs)[1]
	
	## ======================================= ##
	## -- Extracting one treatment
	simp <- simplexs[[i]]	
	smap <- smaps[[i]]
	smapT <- smap[which(smap$theta!=0), ]; smap0 <- smap[which(smap$theta==0), ]
	ar <- ars[[i]]
	
	## ======================================= ##	
	## -- Merge result
	
	if( all(simp[, c(1,2,4)]== smapT[,c(1,2,3)] & smapT[,c(1,2,3)]==ar[,c(1,2,3)] ) ){
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- ASV property
		base <- smapT[,1:5]
		treat = unique(base[,1:2])
		
		summ <- c()
		for(l in 1:nrow(treat)){ #l=1
			
			sub <- Subset(base, treat, l)
			summ <- rbind(summ , sub)
		}
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		stat.summ <- cbind(summ, simplex= simp[5:7],  smapT=smapT[,6:8], smap0=smap0[,6:8], ar[,-c(1,2,3)])			
	}
	
	rho <- gather(stat.summ, method, rho, grep('rho', colnames(stat.summ)))
	rmse <- gather(stat.summ, method, rmse, grep('rmse$', colnames(stat.summ)))
	if(  all(rho[,1:5]==rmse[,1:5], na.rm=TRUE) ){ rr <- cbind(rho, rmse=rmse[,'rmse']) }
	rr[,'method'] <- gsub('\\.rho', '', rr[,'method'])

	## ======================================= ##
	## -- Log scaled RMSE
	stat.summ[,grep('rmse', colnames(stat.summ))] <- apply(stat.summ[,grep('rmse', colnames(stat.summ))]+10e-10, 2, log)
	
	statall <- rbind(statall, cbind(treat=i, stat.summ[,-grep('rmse.weight', colnames(stat.summ))]))
	
	## ======================================= ##
	
}

saveRDS(statall, 'Table/02_population_dynamics_prediction.rds')
############################################################################
## -- Visualizing 

statall <- readRDS('Table/02_population_dynamics_prediction.rds')

pred.fill=c(simplex='royalblue2', smapT='steelblue1', smap0='brown1',  null='darkolivegreen3')
pred.color=c(simplex='royalblue4', smapT='skyblue4', smap0='firebrick3',  null='darkolivegreen4')

## ======================================= ##
## -- Accuracy
 v <- factor(statall$treat, 
			levels=c( 'Soil/medium A', 'Water/medium A',
					  'Soil/medium B', 'Water/medium B', 
					  'Soil/medium C', 'Water/medium C'))		

pdf(sprintf('%s/x=tp_y=accuracy.pdf', dir$figdir), w=15, h=10)
for(k in c('rho', 'rmse$')){#k='rmse$'
	
	
	lf <- gather(statall, key, value, c(grep(sprintf('\\.%s', k), colnames(statall))))
	lf2 <- cbind(lf, group=do.call(rbind, strsplit(lf$key, '\\.')))
	lf2$group.1 <- factor(lf2$group.1, level=names(pred.fill))
	lf2$value[lf2$value>10 | lf2$value == -Inf] <- NA

	pvalsub <- pvalall[ which(pvalall$group.2==k),]
	pvalsub$method <- factor(pvalsub$method, level=names(pred.fill))

	box <- function(x) {
    	v <- c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x))
    	names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    	v
  	}
	
	g <- ggplot(na.omit(lf2))+
		 geom_point(aes(x=as.factor(tp), y=as.numeric(value), group=group.1, color= group.1, fill= group.1), shape=21, 
		 		     position= position_jitterdodge(dodge.width=.8) , size=0.5, alpha=0.8 )+	
		 stat_summary( fun.data=box, geom = "boxplot",width=0.7,size=1, alpha=0.6,
	            	  aes(x=as.factor(tp), y=value, fill= group.1, color= group.1),position=position_dodge(.8) )+
	     geom_text(data= pvalsub, aes(x=as.factor(tp), y=y, label=sig, group=method),position=position_dodge(.8), size=2) +	  
		 facet_wrap( ~ treat, scales='free', nrow=3, dir='h')+
		 scale_fill_manual(values= pred.fill)+
		 scale_color_manual(values=pred.color)+
		 theme_minimal()+
		 theme(panel.border=element_rect(fill=NA, size=1))+
		 labs(x='Forecast time ahead', y=ifelse(k=='rho', 'Rho', 'RMSE (log scale)'), legend='Method')
	plot(g)	 	
}

dev.off()
## ============================================ ##
## -- Dynamics
rm=c('AR')
rm2=c('null')

for(i in names(simplexpc)){#i= names(simplexpc)[3]
	
	## ======================================== ##
    ts <- dlist[[i]]
    stattreat <- statall[statall$treat==i, ]
	
	tmp <- lapply(1:8, function(x){
        y <- apply(ts[libmat[x,1]:libmat[x,2],], 2, scale_2)
        rbind(y, NA)
    })
    tsscale <- do.call(rbind, tmp)
    ## ======================================== ##

	simp <- simplexpc[[i]]
	smapT <- smapTpc[[i]]
	smap0 <- smap0pc[[i]]
	mean <- meanpc[[i]]
	
	ploltlists <- list()
	## ======================================== ##
	for(l in 1:4){ #l=2
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Scatter
	    
	    pred.fill2=c(Simplex='royalblue2', SmapT='steelblue1', Smap0='brown1')
	    pred.color2=c(Simplex='royalblue4', SmapT='skyblue4', Smap0='firebrick3')
	    
		lf <- gather( as.data.frame(cbind(inedx=1:880, replicate=c(rep(1,110),rep(2,110),rep(3,110),rep(4,110),rep(5,110),rep(6,110),rep(7,110),rep(8,110)),
				    Actual=na.omit(tsscale[,l]), tp=simp[, 1],
				    Simplex=simp[, c(l+1)], 'SmapT'=smapT[,c( l+1)],'Smap0'=smap0[,c( l+1)])), key, 'Predicted value', -c(1:4))
		lf <- lf[lf$tp%in%c(1,7), ]
		lf$key <- factor(lf$key, levels=c('Simplex', 'SmapT','Smap0'),
							labels=c('Simplex', expression(paste('S-map (optimized ', italic(theta), ')', sep='')),
										expression(paste('S-map (', italic(theta), '=0)', sep=''))))
		lf$tp <- paste(lf$tp, '-day-ahead', sep='')
		g1 <- ggplot(lf)+
			geom_abline(linetype=2, color='grey50')+
			geom_point(aes(y=Actual, x=`Predicted value`, fill=as.factor(replicate)), shape=23, show.legend=FALSE)+
			facet_grid(key~tp, labeller=label_parsed)+
			scale_fill_manual(values=color$replicate)+
			theme_minimal(base_size=7)+
			theme(panel.border =element_rect(fill=NA, size=0.8))+
			labs(title=colnames(tsscale)[l], y='Observed population size', x='Predicted population size')
			
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Boxplot
		statsub <- stattreat[ which(stattreat $asv==colnames(tsscale)[l]), ]
		statsub <- statsub[,-grep(rm, colnames(statsub))]
		statsub <- statsub[,-grep(rm2, colnames(statsub))]
		
		lf2 <- gather(statsub, key, value, grep('rho', colnames(statsub)))
		lf22 <- cbind(lf2, group=do.call(rbind, strsplit(lf2$key, '\\.')))
		lf22$group.1 <- gsub('s', 'S', lf22$group.1 )
		lf22$group.1 <- factor(lf22$group.1, levels=c('Simplex', 'SmapT','Smap0'),
							labels=c('Simplex', expression(paste('S-map (optimized ', italic(theta), ')', sep='')),
										expression(paste('S-map (', italic(theta), '=0)', sep=''))))
		if(all(is.na(lf22$value))) lf22$value <- 0
		names(pred.color2) <- NULL
		g2 <- ggplot(lf22)+
			  geom_boxplot(aes(x=as.factor(tp), y=value, color=group.1), show.legend=FALSE, alpha=0.7)+
			  facet_wrap(~group.1, ncol=1, labeller=label_parsed)+
			 scale_color_manual(values=pred.color2)+
			 theme_minimal(base_size = 7)+
			 theme(panel.border=element_rect(fill=NA, size=1),
			 	   strip.placement='inside')+
			 labs(x='Time steps in forecasting ', y='Correlation between predicted/observed population size', legend='Method', subtitle='')+
		     scale_y_continuous(breaks=c(-0.5, 0, 0.5), limits=c(-1, 1))

		glist1 <- plot_grid(g1, g2,  nrow=1, align='v', rel_widths=c(0.5, 0.25))
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Dynamics
		lf3 <- gather( data.frame(inedx=rep(1:110, 8), replicate=c(rep('Replicate\ncommunity 1',110),rep('Replicate\ncommunity 2',110),
		                                                           rep('Replicate\ncommunity 3',110),rep('Replicate\ncommunity 4',110),rep('Replicate\ncommunity 5',110),
		                                                           rep('Replicate\ncommunity 6',110),rep('Replicate\ncommunity  7',110),rep('Replicate\ncommunity  8',110)),
				     Actual=na.omit(tsscale[,l]), 
				     Simplex=simp[simp[,1]==1, c(l+1)],'SmapT'=smapT[smapT[,1]==1,c( l+1)],'Smap0'=smap0[smap0[,1]==1,c( l+1)]), key, value, -c(1:2))
		
		g4 <- ggplot()+
			  #geom_rect(aes(xmin=c(1,221,441,661), xmax=c(110,330,550,770), ymin=-Inf, ymax=Inf), fill='snow2', alpha=0.5)+
		      #geom_text( aes(x=c(55, 165, 275, 385, 495, 605, 715, 825), y=max(lf3$value, na.rm=TRUE)-1), label=paste('Replicate', 1:8))+
			  scale_x_continuous(expand=c(0,0))+ 
			  geom_line(data= lf3 ,aes(x= inedx, y=value, group=key, color=key), show.legend=FALSE)+
			  scale_color_manual(values= c(Actual='grey30',pred.fill2))+
			  theme_minimal(base_size = 7)+
			  theme(panel.border=element_rect(fill=NA, size=1))+ 
			  labs(x='Time index', y='Normarized population size', subtitle='1-day-ahead forecasting')+
		    facet_wrap(~replicate, scales='free_x', nrow=2)
		
		ploltlists[[l]] <- plot_grid(glist1, g4, ncol=1, rel_heights=c(0.6, 0.4))	  
	}
    ggsave(sprintf("%s/Statics_result_each_asv_%s.pdf", dir$figdir, gsub('/', '_',i)), marrangeGrob(ploltlists, nrow=3, ncol=2),w=16,h=24)
    svglite::svglite(sprintf("%s/Statics_result_each_asv_%s.svgz", dir$figdir, gsub('/', '_',i)), width=18/2.5, height=27/2.5)
    plot(plot_grid(plotlist=ploltlists, nrow=2))
    dev.off()
}	
