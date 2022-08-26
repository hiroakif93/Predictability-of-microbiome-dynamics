############################################################################
####
#### R script for Fujita (2019)
####
#### Prediction by using simplex projection
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('../')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('ggplot2', 'cowplot', 'RColorBrewer', 'tidyr', 'scales','Rtsne', 'vegan', 'ggrepel'))

# -- Create directory to save
dir <- make.dir('02_Predicting_microbiome_dynamics/Community_dynamics')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')
abruptness <- readRDS('Table/01_Abruptness.rds')
sml <- cbind(sml, abruptness[rownames(sml),])

color <- readRDS('Table/color_palette.rds')
taxa <- readRDS('Table/Taxa_list.rds')

simplexlist <- readRDS('Table/02_06_simplexModel.rds')    
smapTlist <- readRDS('Table/02_07_smapTModel.rds')    
smap0list <- readRDS('Table/02_07_smap0Model.rds')    
meanlist <- readRDS('Table/02_08_meanModel.rds')    
arlist <-readRDS('Table/02_08_ARModel.rds')    

############################################################################
pred.fill=c(simplex='royalblue2', smapT='steelblue1', smap0='firebrick1')

statall <- nmdsl <- nmdsstress <- predictcomm <- c()
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
	
	## ======================================= ##
	## -- Error	 of Absolute value
	errorlist <- comms <- c()	
	for(t in 1:10){#t=5
	    
	    samplerows <- rownames(smlsub[as.numeric(smlsub$time) %in% 21:110, ])
		obstmp <- as.matrix(na.omit(obs))[samplerows,]
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		rownames(smapT)
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

		errors <- cbind(simperror, smapTerror, smap0error)
		colnames(errors)	<- paste(c(names(pred.fill)), '_abs_tp_', t, sep='')
		errorlist <- cbind(errorlist, errors)
		
		comms[[t]] <- cbind(simptp, smapTtmp, smap0tmp)
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
	    
	    errors <- cbind(simperror, smapTerror, smap0error)
	    colnames(errors)	<- paste(c(names(pred.fill)), '_rel_tp_', t, sep='')
	    errorlist <- cbind(errorlist, errors)
	    
	}
	
	#pairs(errors)
	stattmp <- cbind(stattmp, as.data.frame(errorlist)[rownames(smlsub),])
	
	## =============================================================== ##
	## -- Abrupt point index
	sml_re <- c()
	for(r in 1:8){ #r=2
	    sub <- stattmp[stattmp$replicate.id==r, ]
	    sub$abrupt_index <- 1:110-order(sub[,'abruptness_rel_tw5_tp1'], decreasing=TRUE)[1]
	    sml_re <- rbind(sml_re, sub)
	}
	stattmp <- sml_re
	
	## =============================================================== ##
	## Visualising community assembly prediction
	
	pdf(sprintf('%s/predicted_community_assembly_%s.pdf', dir$figdir, gsub('/', '_', i)), w=15)
	for(t in 1:10){#t=1

		#sort(colSums(smapT[smapT[,'tp']==t,-1][441:550,], na.rm=TRUE))
		colnames(simp)[-1] <- colnames(smapT)[-1]  <- colnames(smap0)[-1]  <- colnames(obstmp)
		pred.comms <- cbind(stattmp, actual=as.data.frame(obs)[rownames(stattmp),], simplex=simp[simp[,'tp']==t,-1], 
							smapT = smapT[smapT[,'tp']==t,-1], smap0=smap0[smap0[,'tp']==t,-1])
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
				theme(panel.border=element_rect(size=0.8, fill=NA))+
				xlim(21,110)
				
		lfpred2$value[lfpred2$value<0] <- 0
		rel <- 	ggplot(lfpred2)+
    		    geom_bar(aes(x=as.numeric(time), y=value, fill=group.2), color=NA,
    		             stat='identity', position='fill', show.legend=FALSE, width=1)+
    		    geom_vline(data=lfpred2[which(lfpred2$abrupt_index==0), ], aes(xintercept=as.numeric(time)),
    		               color='red')+
    		    geom_point(data=lfpred2[which(lfpred2$abrupt_index==0), ], aes(x=as.numeric(time), y=1.1 ), fill='grey90', shape=25, size=4)+
    		    facet_grid(group.1~replicate.id, scales='free')+
    		    scale_fill_manual(values=color[[1]][[1]])+
    		    theme_minimal(base_size=26)+
    		    theme(
    		          axis.text.x = element_text(vjust=1, margin=margin(-0.05,0,0,0, unit='inch')),
    		          axis.text.y = element_text(hjust=1, margin=margin(0,0,-0.05,0, unit='inch')))+
    		    labs(title=i, subtitle=sprintf('tp=%s',t))+
    		    scale_x_continuous(expand=c(0,0), limits=c(21,110),breaks=seq(30, 90, 30))+
    		    scale_y_continuous(expand=expansion(mult = c(0, .05)),breaks=seq(0.25, 0.75, 0.25))

		plot(abs); plot(rel)		
	}
	dev.off()

	## =============================================================== ##
    ## -- Visualizing by NMDS
	obsdf <- as.matrix(as.data.frame(obs)[samplerows,])
	
	nmdsdf <- c()
	for(tp in c(1,7)){
	    
	    df2 <- as.matrix(rbind(obsdf/rowSums(obsdf), 
	                           smapT[smapT[,'tp']==tp,][samplerows,-1]/rowSums(smapT[smapT[,'tp']==tp,][samplerows,-1], na.rm=TRUE))) 
	    df2[is.na(df2)] <- 0
	    
	    nmds <- metaMDS(df2+0.00001, parallel=20)
	    nmdstmp <- cbind(tp=tp, sample=c(rownames(obsdf), rownames(smapT[smapT[,'tp']==tp,][samplerows,-1])), 
	                    type=c(rep('obs', nrow(nmds$point)/2), rep('pred', nrow(nmds$point)/2)), 
	                    stattmp[rownames(df2),], nmds$point)

	    nmdsstress <- rbind(nmdsstress, data.frame(treat=i, tp=tp, nmds$stress))
	    nmdsdf <- rbind(nmdsdf, nmdstmp)
	}

	## =============================================================== ##
	
	statall <- rbind(statall, stattmp)
	nmdsl <- rbind(nmdsl, nmdsdf)
}

saveRDS(statall, 'Table/02_Error_community_prediction.rds')
saveRDS(predictcomm, 'Table/02_Community_dynamics_prediction.rds')
saveRDS(nmdsl, 'Table/02_Prediction_NMDS.rds')

## ======================================= ##

data <- readRDS('Table/02_Prediction_NMDS.rds')

glist <- c()
for(i in names(dlist)[-7]){#i=names(dlist)[2]
    
    sub <- data[data$treat1==i, ]
    
    nmdsg1 <- ggplot(sub)+
                geom_line(aes(x=MDS1, y=MDS2, group=sample), color='grey50')+
                geom_point(aes(x=MDS1, y=MDS2, shape=type, fill=as.factor(replicate.id) ), color='grey40', size=2)+
                scale_fill_manual(values=color$replicate)+
                scale_shape_manual(values=c(21,24))+
                theme_minimal()+
                theme(panel.border = element_rect(size=1, fill=NA),
                      legend.background = element_blank(),
                      legend.text = element_text(hjust=0))+ 
                guides(fill = guide_colourbar(title.position = "top"),
                       shape = guide_legend(title.position = "top"))+
                coord_fixed(ratio=diff(range(sub$MDS1))/diff(range(sub$MDS2)))+
                labs(title=i, x='Axis 1', y='Axis 2', fill='Time step\nfrom maximal change')+
                facet_wrap(~tp)
                
    glist[[i]] <- nmdsg1
}

ggsave(plot=plot_grid(plotlist=glist), 
       filename=sprintf('%s/nmds_relative.pdf', dir$figdir), w=23, h=10)
       

