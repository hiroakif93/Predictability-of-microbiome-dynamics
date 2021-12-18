############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot', 'RColorBrewer', 'doParallel', 'Rtsne', 'ggsnippets', 'cluster', 'fpc'))

# -- Create directory to save
dir <- make.dir('03_dynamics_indices/02_assembly')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/03_matrixList.rds')
sml <- readRDS('Table/03_sample_info.rds')

taxa <- readRDS('Table/Taxa_list_02.rds')
color <- readRDS('Table/03_color_palette_v2.rds')

############################################################################
nmdsl <- statall <- tsnel <- c()

for(i in names(dlist)[-7]){ #i=names(dlist)[5]
	
	## ============================================================ ##
	## -- Extracting one treatment matrix
	
	smlsub <- sml[sml$treat1==i, ]
	ts <- dlist[[i]]
	rel.ts <- ts/rowSums(ts)
	
	## ======================================= ##
	## -- Clustering 
	smltmp <- c()
	for(r in 1:8){ #r=1
		
		smlrep <- smlsub[smlsub$replicate.id==r,]
		
		tssub <- ts[rownames(smlrep), ]
		
		smlrep$cluster <- cutree(hclust( vegdist(tssub/rowSums(tssub)) ), k=2)
		
		smltmp <- rbind(smltmp, smlrep)
	}
	
	## ======================================= ##
	## -- Visualizing community assembly
	pdf(sprintf('%s/community_assembly_%s.pdf', dir$figdir, gsub('/', '_', i)), w=15)
	for(l in names(color$asv)){ #l=colnames(taxa)[3]
		
		agg.ts <- t(Taxa.mat(ts[rownames(smlsub),], taxa, l))
		
		lf <- gather(cbind(smlsub, agg.ts), 
					 key, value, -c(1:(ncol(smlsub))))
		
		lf$key <- factor(lf$key, levels=names(sort(colSums(agg.ts))) )
		gtmp <- ggplot(lf)+
			  geom_area(aes(x=as.numeric(time), y=value, fill=key), stat='identity')+
			  scale_fill_manual( values= color$asv[[l]][colnames(agg.ts)] )
		gl <- g_legend(gtmp)	  	

  
		g1 <- ggplot(lf)+
			  geom_area(aes(x=as.numeric(time), y=value, fill=key), color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
			  facet_wrap(~replicate.id,scales='free')+
			  scale_fill_manual( values= color$asv[[l]][colnames(agg.ts)] ) + 
	          theme_minimal(base_size=15)+
	          labs(x= "Day", title=i, subtitle=gsub('ID', 'ASV', l), y= "Number of DNA coies")+
	          guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		relts <- agg.ts/rowSums(agg.ts)
		lf <- gather(cbind(smlsub, agg.ts), 
					 key, value, -c(1:(ncol(smlsub))))
		
		lf$key <- factor(lf$key, levels=names(sort(colSums(relts))) )
  
		g2 <- ggplot(lf)+
			  geom_area(aes(x=as.numeric(time), y=value, fill=key), color='grey30', 
			  			stat='identity', show.legend=FALSE, size=0.1, position='fill')+
			  facet_wrap(~replicate.id,scales='free')+
			  scale_fill_manual( values= color$asv[[l]][colnames(agg.ts)] ) + 
	          theme_minimal(base_size=15)+
	          labs(x= "Day", title=i, subtitle=gsub('ID', 'ASV', l), y= "Number of DNA coies")+
	          guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))
		
		plot(g1)		
		plot(g2)
		plot(gl)
		
	}
	dev.off()
	
	## ======================================= ##
	## -- NMDS
	abs <- ts[rownames(smlsub),]
	rel <- abs/rowSums(abs)
	
	absnmds <- metaMDS(abs, distance='bray', k=2, trymax=20, parallel=detectCores())
	relnmds <- metaMDS(rel, distance='bray', k=2, trymax=20, parallel=detectCores())
	outlier <- c()
	if( any( apply(absnmds $points, 2, sd)>1 ) ){
		
		while( any( apply(absnmds$points, 2, sd)>1 )){	
			out <- which(absnmds $points[, which(apply(absnmds $points, 2, sd)>1)] == max( absnmds $points[, which(apply(absnmds $points, 2, sd)>1)] ))
			outlier <- c(outlier, rownames(ts)[out])
			absnmds <- metaMDS(abs[-out, ], distance='bray', k=2, trymax=20, parallel=detectCores())
			relnmds <- metaMDS(rel[-out, ], distance='bray', k=2, trymax=20, parallel=detectCores())
		}
	}
	
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	## -- tSNE
	abstsne <- Rtsne(abs, check_duplicates = FALSE)
	reltsne <- Rtsne(rel, check_duplicates = FALSE)
	
	## ======================================= ##
									
	lf2 <- rbind(cbind(smlsub, type='Absolute', absnmds$points, tSNE= abstsne$Y), 
				 cbind(smlsub, type='Relative', relnmds$points, tSNE= reltsne$Y))
	
	gnmds <- ggplot(lf2)+
			geom_point(aes(x=MDS1, y=MDS2, color=as.numeric(time), shape=as.factor(replicate.id)))+
			facet_wrap(~type, scales='free')+
			scale_color_gradientn(colors=brewer.pal(11,'Spectral')[-6])+
			theme_minimal(base_size=15)+
			theme(panel.border=element_rect(fill=NA, size=0.8))+
			scale_shape_manual(values=c(16:18,15, 7:10))+
			labs(x='Axis 1', y='Axis 2', caption=c(sprintf('Absolute stress : %s\nRelative stress : %s\n', round(absnmds$stress, digits=3), round(relnmds$stress, digits=3) )),
				 title=i)
	
	nmdsl[[i]] <- gnmds		
	
	gtsn1 <- ggplot(lf2)+
			 geom_point(aes(x=tSNE.1, y=tSNE.2, color=as.numeric(time), shape=as.factor(replicate.id)))+
			 facet_wrap(~type, scales='free')+
			 scale_color_gradientn(colors=brewer.pal(11,'Spectral')[-6])+
			 theme_minimal(base_size=15)+
			 theme(panel.border =element_rect(fill=NA, size=0.8))+
			 scale_shape_manual(values=c(16:18,15, 7:10))+
			 labs(x='Axis 1', y='Axis 2', title=i)
	gtsn2 <- ggplot(lf2)+
			 geom_point(aes(x=tSNE.1, y=tSNE.2, color=as.factor(replicate.id), shape=as.factor(replicate.id)))+
			 facet_wrap(~type)+
			 scale_color_manual(values=color$replicate)+
			 theme_minimal()+
			 theme(panel.border =element_rect(fill=NA, size=0.8))+
			 scale_shape_manual(values=c(16:18,15, 7:10))
			 
	tsnel[[i]] <- gtsn1
	## ======================================= ##
		
	statall <- rbind(statall, cbind(smltmp, abs=absnmds$points[rownames(smlsub),], rel=relnmds$points[rownames(smlsub),],
									abs.tsne= abstsne$Y, rel.tsne = reltsne$Y ) )
	
}
	 

saveRDS(statall,'Table/smaple_info_03_02.rds')

pdf(sprintf('%s/nmds.pdf', dir$figdir), w=17)
for(i in 1:length(nmdsl)){ plot(nmdsl[[i]])}
dev.off()

pdf(sprintf('%s/tsne.pdf', dir$figdir), w=15)
for(i in 1:length(tsnel)){ plot(tsnel[[i]])}
dev.off()
############################################################################
