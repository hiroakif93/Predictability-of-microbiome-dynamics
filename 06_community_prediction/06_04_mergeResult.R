############################################################################
####
#### R script for Fujita (2019)
####
#### Prediction by using simplex projection
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('ggplot2', 'cowplot', 'RColorBrewer', 'tidyr', 'scales','Rtsne', 'vegan', 'ggrepel'))

# -- Create directory to save
dir <- make.dir('06_community_prediction/06_04_output')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/03_matrixList.rds')
libmat=dlist[[7]]
sml  <- readRDS('Table/03_06_sample_info.rds')
color <- readRDS('Table/03_color_palette.rds')
taxa <- readRDS('Table/Taxa_list_02.rds')

simplexlist <- readRDS('Table/06_01_simplex_predicted_community')    
smapTlist <- readRDS('Table/06_02_smapT_predicted_community')    
smap0list <- readRDS('Table/06_02_smap0_predicted_community')    
meanlist <- readRDS('Table/06_03_mean_predicted_community')    
arlist <-readRDS('Table/06_03_AR_predicted_community')    

# -- Decided regularized parameter range
E.range <- 1:20
time.horizon <- 1:10
tp=c(1:10)

############################################################################
pred.fill=c(simplex='royalblue2', smapT='steelblue1', smap0='firebrick1', AR='orange1', null='darkolivegreen3')

statall <- nmdsl <- tsnel <- predictcomm <- absl <- rell <- nmdsdf1 <- nmdsdf2 <- c()
for(i in names(simplexlist)){#i=names(simplexlist)[2]
	
	## ======================================= ##
	## -- Actual / Predicted communities
	obs <- dlist[[i]]
	smlsub <- sml[sml$treat1==i, ]

	simp <- simplexlist[[i]]
	
	smapT <- smapTlist[[i]]
	smapT[smapT <0] <- 0
	
	smap0 <- smap0list[[i]]
	smap0[smap0 <0] <- 0
	
	ar <- arlist[[i]]
	ar[ar <0] <- 0
	
	mean <- meanlist[[i]]
	
	## ======================================= ##
	## -- Error	 of Absolute value
	errorlist <- comms <- c()	
	for(t in 1:10){#t=5
	    
	    samplerows <- rownames(smlsub[as.numeric(smlsub$time) %in% 21:110, ])
		obstmp <- as.matrix(na.omit(obs))[samplerows,]
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
		simptp <- simp[which(simp[,1]==t),-1][samplerows,]
		simperror <- rowSums(abs(obstmp-simptp)) / ( rowSums(obstmp) + rowSums(simptp) )
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
				
		smapTtmp <- as.matrix(smapT[which(smapT[,1]==t),-1])[samplerows,]		
		smapTtmp[is.na(smapTtmp)] <- 0
		smapTerror <- rowSums(abs(obstmp- smapTtmp)) / ( rowSums(obstmp) + rowSums(smapTtmp) )
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
				
		smap0tmp <- as.matrix(smap0[which(smap0[,1]==t),-1])[samplerows,]		
		smap0tmp[is.na(smap0tmp)] <- 0
		smap0error <- rowSums(abs(obstmp- smap0tmp)) / ( rowSums(obstmp) + rowSums(smap0tmp) )
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		artmp <- ar[which(ar[,1]==t),-1][samplerows,]
		artmp[is.na(artmp)] <- 0
		ARerror <-rowSums(abs(obstmp-artmp)) / ( rowSums(obstmp) + rowSums(artmp) )
		meanerror <- rowSums(abs(obstmp-mean[which(mean[,1]==t),-1][samplerows,])) / ( rowSums(obstmp) + rowSums(mean[which(mean[,1]==t),-1][samplerows,]) )		

		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

		errors <- cbind(simperror, smapTerror, smap0error, ARerror, meanerror, meanerror-simperror, delta_error=meanerror-smapTerror)
		colnames(errors)	<- paste(c(names(pred.fill), 'delta_simplex_error', 'delta_smap_error'), '_abs_tp_', t, sep='')
		errorlist <- cbind(errorlist, errors)
		
		comms[[t]] <- cbind(simptp, smapTtmp, smap0tmp, artmp, mean[which(mean[,1]==t),-1][samplerows,])
	}
	predictcomm[[i]] <- comms
	stattmp <- cbind(smlsub, index=1:nrow(smlsub), as.data.frame(errorlist)[rownames(smlsub),])
	
	## ============================================================================== ##
	## -- Error	 of Relative value
	errorlist <- c()	
	for(t in 1:10){#t=5
	    
	    samplerows <- rownames(smlsub[as.numeric(smlsub$time) %in% 21:110, ])
	    obstmp <- as.matrix(na.omit(obs))[samplerows,]
	    obstmp <- obstmp/rowSums(obstmp)
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    
	    simptp <- simp[which(simp[,1]==t),-1][samplerows,]
	    simptp <- simptp/rowSums(simptp)
	    simperror <- rowSums(abs(obstmp-simptp)) / ( rowSums(obstmp) + rowSums(simptp) )
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    
	    smapTtmp <- as.matrix(smapT[which(smapT[,1]==t),-1])[samplerows,]		
	    smapTtmp[is.na(smapTtmp)] <- 0
	    smapTtmp <- smapTtmp/rowSums(smapTtmp)
	    smapTerror <- rowSums(abs(obstmp- smapTtmp)) / ( rowSums(obstmp) + rowSums(smapTtmp) )
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    
	    smap0tmp <- as.matrix(smap0[which(smap0[,1]==t),-1])[samplerows,]		
	    smap0tmp[is.na(smap0tmp)] <- 0
	    smap0tmp <- smap0tmp/rowSums(smap0tmp)
	    smap0error <- rowSums(abs(obstmp- smap0tmp)) / ( rowSums(obstmp) + rowSums(smap0tmp) )
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    artmp <- ar[which(ar[,1]==t),-1][samplerows,]
	    artmp[is.na(artmp)] <- 0
	    artmp <- artmp/rowSums(artmp)
	    ARerror <-rowSums(abs(obstmp-artmp)) / ( rowSums(obstmp) + rowSums(artmp) )
	    
	    meantmp <- mean[which(mean[,1]==t),-1][samplerows,]
	    meantmp <- meantmp/rowSums(meantmp)
	    meanerror <- rowSums(abs(obstmp-meantmp)) / ( rowSums(obstmp) + rowSums(meantmp) )		
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    
	    errors <- cbind(simperror, smapTerror, smap0error, ARerror, meanerror, meanerror-simperror, delta_error=meanerror-smapTerror)
	    colnames(errors)	<- paste(c(names(pred.fill), 'delta_simplex_error', 'delta_smap_error'), '_rel_tp_', t, sep='')
	    errorlist <- cbind(errorlist, errors)
	    
	}
	
	#pairs(errors)
	stattmp <- cbind(stattmp, as.data.frame(errorlist)[rownames(smlsub),])
	## =============================================================== ##
	## Visualising community assembly prediction
	
	pdf(sprintf('%s/pred_community_%s.pdf', dir$figdir, gsub('/', '_', i)), w=15)
	for(t in 1:10){#t=1

		sort(colSums(smapT[smapT[,'tp']==t,-1][441:550,], na.rm=TRUE))
		colnames(simp)[-1] <- colnames(smapT)[-1]  <- colnames(ar)[-1]  <- colnames(smap0)[-1]  <- colnames(obstmp)
		pred.comms <- cbind(stattmp, actual=as.data.frame(obs)[rownames(stattmp),], simplex=simp[simp[,'tp']==t,-1], 
							smapT = smapT[smapT[,'tp']==t,-1], AR=ar[ar[,'tp']==t,-1], smap0=smap0[smap0[,'tp']==t,-1])
		lfpred <- gather(pred.comms, key, value, -c(1:ncol(stattmp)))
		lfpred2 <- cbind(lfpred, group=do.call(rbind, strsplit(lfpred$key, '\\.')))
		lfpred2$group.1 <- factor(lfpred2$group.1, levels=c('actual',names(pred.fill)))
		abs <- ggplot(lfpred2)+
			geom_bar(aes(x=as.numeric(time), y=value, fill=group.2), stat='identity', show.legend=FALSE)+
			#geom_vline(xintercept=c(errorindex[which(errorindex$index==0), 'time']))+
			facet_grid(group.1~replicate.id, scales='free')+
			scale_fill_manual(values=color$asv[[1]])+
			theme_minimal()+
			theme(panel.border=element_rect(size=0.8, fill=NA))+
			xlim(21,110)+
			labs(title=i, subtitle=sprintf('tp=%s',t))
				
		lfpred2$value[lfpred2$value<0] <- 0
		rel <- 	ggplot(lfpred2)+
				geom_bar(aes(x=as.numeric(time), y=value, fill=group.2), stat='identity', position='fill', show.legend=FALSE)+
				#geom_vline(xintercept=c(errorindex[which(errorindex$index==0), 'time']))+
				facet_grid(group.1~replicate.id, scales='free')+
				scale_fill_manual(values=color$asv[[1]])+
				theme_minimal()+
				theme(panel.border=element_rect(size=0.8, fill=NA))+
				xlim(21,110)+
				labs(title=i, subtitle=sprintf('tp=%s',t))
		plot(abs); plot(rel)		
	}
	dev.off()
	
	pdf(sprintf('%s/pred_community_%s_ver2.pdf', dir$figdir, gsub('/', '_', i)), w=15)
	for(t in 1:10){#t=1

		#sort(colSums(smapT[smapT[,'tp']==t,-1][441:550,], na.rm=TRUE))
		colnames(simp)[-1] <- colnames(smapT)[-1]  <- colnames(ar)[-1]  <- colnames(smap0)[-1]  <- colnames(obstmp)
		pred.comms <- cbind(stattmp, actual=as.data.frame(obs)[rownames(stattmp),], simplex=simp[simp[,'tp']==t,-1], 
							smapT = smapT[smapT[,'tp']==t,-1], smap0=smap0[smap0[,'tp']==t,-1], Null=mean[mean[,'tp']==t,-1])
		lfpred <- gather(pred.comms, key, value, -c(1:ncol(stattmp)))
		lfpred2 <- cbind(lfpred, group=do.call(rbind, strsplit(lfpred$key, '\\.')))
		lfpred2$group.1 <- factor(lfpred2$group.1, levels=c('actual',names(pred.fill)))
		
		abs <- ggplot(lfpred2)+
				geom_bar(aes(x=as.numeric(time), y=value, fill=group.2), stat='identity', show.legend=FALSE, width=1)+
				geom_vline(xintercept=c(lfpred2[which(lfpred2$abruptness_rel_tw5_tp1==max(lfpred2$abruptness_rel_tw5_tp1, na.rm=TRUE)), 'time']),
						  color='red')+
				facet_grid(group.1~replicate.id, scales='free')+
				scale_fill_manual(values=color$asv[[1]])+
				theme_minimal()+
				theme(panel.border=element_rect(size=0.8, fill=NA),
				      text=element_text(family='Arial'))+
				xlim(21,110)
				
		lfpred2$value[lfpred2$value<0] <- 0
		rel <- 	ggplot(lfpred2)+
				geom_bar(aes(x=as.numeric(time), y=value, fill=group.2), stat='identity', position='fill', show.legend=FALSE, width=1)+
				geom_vline(xintercept=c(lfpred2[which(lfpred2$abruptness_rel_tw5_tp1==max(lfpred2$abruptness_rel_tw5_tp1, na.rm=TRUE)), 'time']),
						  color='red')+
				facet_grid(group.1~replicate.id, scales='free')+
				scale_fill_manual(values=color$asv[[1]])+
				theme_minimal()+
				theme(panel.border=element_rect(size=0.8, fill=NA))+
				labs(title=i, subtitle=sprintf('tp=%s',t))+
				scale_x_continuous(expand=c(0,0), limits=c(21,110))+scale_y_continuous(expand=c(0,0))
		plot(abs); plot(rel)		
	}
	dev.off()
	
	absl[[i]] <- abs
	rell[[i]] <- rel	
	
	## =============================================================== ##
    ## -- Visualizing by NMDS and tSNE
	
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	stattmp$abruptIndex <- NA
	for(r in 1:8){#r=1
		
		reprow <- which(stattmp[,'replicate.id']==r)
		maxabrupt <- which(stattmp[reprow, 'abruptness_rel_tw5_tp1']==max(stattmp[reprow, 'abruptness_rel_tw5_tp1'], na.rm=TRUE))
		stattmp[reprow, 'abruptIndex'] <- c(1:110)-maxabrupt
		
	}
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	## -- Absolute
	obsdf <- as.matrix(as.data.frame(obs)[samplerows,])
	df <- as.matrix(rbind(obsdf, smapT[samplerows,-1])); df[is.na(df)] <- 0
	
	nmds <- metaMDS(df+0.00001, parallel=7)
	tsne <- Rtsne(df, check_duplicates=FALSE)$Y
	rownames(tsne) <- rownames(df)
	nmds.sample <- cbind(sample=c(rownames(obsdf[rownames(df),]), rownames(obsdf[rownames(df),])), 
	                      type=c(rep('obs', nrow(nmds$point)/2), rep('pred', nrow(nmds$point)/2)), 
	                      rbind(stattmp[rownames(df),], stattmp[rownames(df),] ), nmds$point, tsne=tsne)
	
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	## -- Relative
	df2 <- as.matrix(rbind(obsdf/rowSums(obsdf), smapT[smapT[,'tp']==1,][samplerows,-1]/rowSums(smapT[smapT[,'tp']==1,][samplerows,-1], na.rm=TRUE))); 
	df2[is.na(df2)] <- 0
	nmds2 <- metaMDS(df2+0.00001, parallel=7)

	tsne2 <- Rtsne(df2, check_duplicates=FALSE)$Y
	rownames(tsne2) <- rownames(df2)

	nmds.sample2 <- cbind(sample=c(rownames(obsdf[rownames(df2),]), rownames(obsdf[rownames(df2),])), 
						type=c(rep('obs', nrow(nmds2$point)/2), rep('pred', nrow(nmds2$point)/2)), 
						rbind(stattmp[rownames(df2),], stattmp[rownames(df2),] ), nmds2$point, tsne=tsne2)

	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	
	deltamin <- min(nmds.sample$delta_smap_error_abs_tp_1, na.rm=TRUE)
	deltamax <- max(nmds.sample$delta_smap_error_abs_tp_1, na.rm=TRUE)
	nmdsg1 <- ggplot(nmds.sample)+
			  geom_line(aes(x=MDS1, y=MDS2, group=sample), color='grey50')+
			  geom_point(aes(x=MDS1, y=MDS2, shape=type, fill=as.numeric(time)), color='grey40', size=3)+
	          scale_fill_gradientn(colors=brewer.pal(11, 'RdBu')[-c(5,6)], na.value=NA, 
	                         values=rescale(c(deltamin, 0, deltamax)))+
	          scale_shape_manual(values=c(21,24))+
	          theme_bw()+
	          theme(text=element_text(family='Arial'))+
	          coord_fixed(ratio=diff(range(nmds.sample$MDS1))/diff(range(nmds.sample$MDS2)))+
	          labs(title=i, x='Axis 1', y='Axis 2')

	tsneg1 <- ggplot(nmds.sample)+
			  geom_point(aes(x= tsne.1, y=tsne.2, shape=type, fill=delta_smap_error_abs_tp_1), color='grey40')+
			  geom_line(aes(x= tsne.1, y=tsne.2, group=sample), color='grey50')+
    	      scale_fill_gradientn(colors=brewer.pal(11, 'RdBu')[-c(5,6)], na.value=NA, 
    	                         values=rescale(c(deltamin, 0, deltamax)))+
    	      scale_shape_manual(values=c(21,24))+
    	      theme_bw()+
    	      coord_fixed(ratio=diff(range(nmds.sample$tsne.1))/diff(range(nmds.sample$tsne.2)))+
	    labs(title=i)
	
	nmdsg2 <- ggplot(nmds.sample2)+
			  geom_line(aes(x=MDS1, y=MDS2, group=sample), color='grey50')+
    	      geom_point(aes(x=MDS1, y=MDS2, shape=type, fill=as.numeric(abruptIndex)), color='grey40', size=3)+
    	      scale_fill_gradientn(colors=brewer.pal(11, 'RdBu')[-c(5,6)], na.value=NA, 
    	                         values=rescale(c(min(nmds.sample2$abruptIndex, na.rm=TRUE), 0, max(nmds.sample2$abruptIndex, na.rm=TRUE))))+
    	      scale_shape_manual(values=c(21,24))+
    	      theme_bw()+
    	      coord_fixed(ratio=diff(range(nmds.sample2$MDS1))/diff(range(nmds.sample2$MDS2)))+
	    	  theme(text=element_text(family='Arial'))+
	          labs(title=i, x='Axis 1', y='Axis 2')

	tsneg2 <- ggplot(nmds.sample2)+
			 geom_point(aes(x= tsne.1, y=tsne.2, shape=type, color=as.numeric(time)))+
	    labs(title=i)
	

	nmdsl[[i]] <- plot_grid(nmdsg1, nmdsg2)
	tsnel[[i]] <- plot_grid(tsneg1, tsneg2)
	
	nmdsdf1 <- rbind(nmdsdf1, nmds.sample)
	nmdsdf2 <- rbind(nmdsdf2, nmds.sample2)
	
	## =============================================================== ##
	
	statall <- rbind(statall, stattmp)
	
}

saveRDS(statall, 'Table/06_04_sample_info.rds')
saveRDS(predictcomm, 'Table/06_04_predict_community.rds')
saveRDS(nmdsdf1, 'Table/06_04_nmds_abs.rds')
saveRDS(nmdsdf2, 'Table/06_04_nmds_rel.rds')

pdf(sprintf('%s/nmds.pdf', dir$figdir), w=17)
for(l in 1:length(nmdsl)){ plot(nmdsl[[l]])}
dev.off()

pdf(sprintf('%s/tsne.pdf', dir$figdir), w=17)
for(l in 1:length(tsnel)){ plot(tsnel[[l]])}
dev.off()

## ======================================= ##

nmdsdf1 <- readRDS('Table/06_04_nmds_abs.rds')
nmdsdf2 <- readRDS('Table/06_04_nmds_rel.rds')

glist <- c()
for(i in names(dlist)[-7]){
    
    sub <- nmdsdf2[nmdsdf2$treat1==i, ]
    
    deltamin <- min(sub$delta_smap_error_abs_tp_1, na.rm=TRUE)
    deltamax <- max(sub$delta_smap_error_abs_tp_1, na.rm=TRUE)
    
    nmdsg1 <- ggplot(sub)+
                geom_line(aes(x=MDS1, y=MDS2, group=sample), color='grey50')+
                geom_point(aes(x=MDS1, y=MDS2, shape=type, fill=as.numeric(abruptIndex)), color='grey40', size=3)+
                scale_fill_gradientn(colors=brewer.pal(11, 'RdBu')[-c(5,6)], na.value=NA, 
                                     values=rescale(c(deltamin, 0, deltamax)))+
                scale_shape_manual(values=c(21,24))+
                theme_minimal(base_size)+
                theme(panel.border = element_rect(size=1, fill=NA),
                        text=element_text(family='Arial'),
                      legend.background = element_blank(),
                      legend.text = element_text(hjust=0))+ 
                guides(fill = guide_colourbar(title.position = "top"),
                       shape = guide_legend(title.position = "top"))+
                coord_fixed(ratio=diff(range(sub$MDS1))/diff(range(sub$MDS2)))+
                labs(title=i, x='Axis 1', y='Axis 2', fill='Time step\nfrom maximal change')
    glist[[i]] <- nmdsg1
}


ggsave(plot=plot_grid(plotlist=glist), 
       filename=sprintf('%s/nmds_relative.tiff', dir$figdir), w=14, h=10)
