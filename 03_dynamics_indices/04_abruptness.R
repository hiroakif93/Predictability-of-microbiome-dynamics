############################################################################
####
#### R script for Fujita (2019)
####
#### Abruptness
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot', 'vegetarian', 'scales'))

# -- Create directory to save
dir <- make.dir('03_dynamics_indices/04_abruptness')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS("Table/03_matrixList.rds")
sml <- readRDS('Table/03_03_sample_info.rds')
color <- readRDS('Table/03_color_palette.rds')
taxa <- readRDS("Table/Taxa_list_02.rds")

time.step = c(1:10)
previous=5
timewindow=c(3,5,7)

for.parallel(8)
############################################################################

abruptness <- all <- c()

pdf(sprintf('%s/abruptness.pdf', dir$figdir), w=15,h=10)
for(i in names(dlist)[-7]){ # i=names(dlist)[2]
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
	
    ts.sub <- na.omit(dlist[[i]]) 
    smlsub <- sml[which(sml $treat1==i),]
    
 	## ========================================================== ##
 	## -- Beta diversity against replicate
 	
 	smlsub$bray_against_replicate <- NA; smlsub$beta.q0_against_replicate <- NA;
	for(n in 1:110){#n=50
		
		smlday <- smlsub[ which(as.numeric(smlsub$time)==n), ]
		
		bray <- colMeans( as.matrix(vegdist(ts.sub[rownames(smlday), ], method='bray')) )
		jost <- matrix(NA, ncol=8, nrow=8, dimnames=list(rownames(smlday), rownames(smlday)))
		for(m in 1:8){
			for(l in 1:8){
				jost[m, l] <- d(ts.sub[c(colnames(jost)[m], colnames(jost)[l]), ], lev='beta')
			}
		}

		smlsub[names(bray), ]$bray_against_replicate <- bray
		smlsub[names(colMeans(jost)), ]$beta.q0_against_replicate <- colMeans(jost)
	}
	
	## ========================================================== ##
    ## -- Bray-curtis dissimilarity between previous and future state 

    abruptness <- c()
	for(tw in timewindow){
		for(tt in time.step){ # tt=5
				
		    tmp.unpred <- foreach(r=1:8, .packages=c('vegan', 'vegetarian')) %dopar%{ #r=4
		        
		        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		        ## -- Extracting one replicaate
		        sample.replicate <- smlsub[ which(smlsub $replicate.id==r), ]
		        
		        ts.abs.sub <- ts.sub[rownames(sample.replicate), ]
		        ts.rel.sub <- ts.sub[rownames(sample.replicate), ]/rowSums(ts.sub[rownames(sample.replicate), ])
	     
		        unp.tmp <- matrix(NA, nrow=110, ncol=1*3)
		        
		        colnames(unp.tmp) <- c(sprintf('abruptness_jost_tw%s_tp%s', tw, tt ),
		                               sprintf('abruptness_rel_tw%s_tp%s', tw,tt ) ,
		                               sprintf('abruptness_abs_tw%s_tp%s',tw, tt ))		        
		        
	        	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		        ## -- Calculation of abruptness
		        
		        for(day in (tw+1):(110-tw-tt)){ #day =1
		            
		            future.id <- rownames(sample.replicate)[which(as.numeric(sample.replicate$time) %in% (day+tt): (day+ tw))]
		            previous.id <- rownames(sample.replicate)[which(as.numeric(sample.replicate$time) %in% (day-1): (day-tw))]
		            
		            ## ========================================================== ##
					## -- Absolute value community
					future.ab <- colMeans( ts.abs.sub[future.id, ] )
					previous.mean.ab <- colMeans(ts.abs.sub[previous.id, ])
					
		            ## -- Bray-Curtis dissimilarity between abundance at t+1 and mean abundance t-time.window+1 : t
		            jost <- d(rbind(future.ab, previous.mean.ab), lev='beta')
		            
		            ## -- Bray-Curtis dissimilarity between abundance at t+1 and mean abundance t-time.window+1 : t
		            abs <- vegdist( rbind(previous.mean.ab, future.ab), method='bray', na.rm=TRUE )
		            
		            ## ========================================================== ##
					## -- Relative abundance community
					future.ab <- colMeans(ts.rel.sub[future.id, ])
					previous.mean.ab <- colMeans(ts.rel.sub[previous.id, ])
					
		            ## -- Bray-Curtis dissimilarity between abundance at t+1 and mean abundance t-time.window+1 : t
		            relative <- vegdist( rbind(previous.mean.ab, future.ab), method='bray', na.rm=TRUE )
		            
		            unp.tmp[day, ] <- c(jost, relative, abs)
		        }
	        	
	        	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		        
		        rownames(unp.tmp) <- rownames(sample.replicate) 
		        unp.tmp
		        
		    } 
		    
		    abruptness <- cbind(abruptness, do.call(rbind, tmp.unpred))   
	    }
	}    
		
	## ========================================================== ##

    smlmerge <- cbind(smlsub, abruptness)
	all <- rbind(all, smlmerge)
	
    ## ========================================================== ##
    lfarea <- gather( cbind(smlmerge, #tmp.domi,
                            ts.sub[rownames(smlmerge),]), key, value, -c(1:(ncol(smlmerge)))) 
    lfarea$key <- factor( lfarea$key, levels= names(sort(colSums(ts.sub))) )
    sec.multip <- max(abruptness[rownames(smlsub),'abruptness_rel_tw5_tp1'], na.rm=TRUE)/max(ts.sub, na.rm=TRUE)
    ggarea <- ggplot(lfarea)+
	            geom_area( aes(x=as.numeric(time), y=value, group=key, fill=key ))+
	            geom_point( aes(x=as.numeric(time), y= rescale(abruptness_rel_tw5_tp1, to=c(0, max(ts.sub))), group=key), 
	                        shape=23, fill='red', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(abruptness_rel_tw5_tp1, to=c(0, max(ts.sub))), group=key), 
	                   color='red', alpha=0.2, size=0.2 ) +
	            geom_point( aes(x=as.numeric(time), y= rescale(abruptness_abs_tw5_tp1, to=c(0, max(ts.sub))), group=key), 
	                        shape=23, fill='blue', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(abruptness_abs_tw5_tp1, to=c(0, max(ts.sub))), group=key), 
	                   color='blue', alpha=0.2, size=0.2 ) + 
	            facet_wrap( ~ replicate.id, scales='free') +
	            scale_fill_manual( values= color$asv[[1]] ) +
	            labs(x="Day", title=i, y= "Number of DNA coies",
	                 subtitle='red=relative, blue=absolute')+
	            guides(fill=guide_legend(title='ASV', nrow=2, reverse=TRUE))+ 
	            theme_minimal()+
	            theme(legend.position='') +
	            scale_y_continuous(sec.axis=sec_axis(~.* sec.multip, name='Abruptness (time window 5, tp 1)')) 
    plot(ggarea)
    ggarea <- ggplot(lfarea)+
	            geom_area( aes(x=as.numeric(time), y=value, group=key, fill=key ), position='fill')+
	            geom_point( aes(x=as.numeric(time), y= rescale(abruptness_rel_tw5_tp1, to=c(0, 1)), group=key), 
	                        shape=23, fill='red', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(abruptness_rel_tw5_tp1, to=c(0, 1)), group=key), 
	                   color='red', alpha=0.2, size=0.2 ) +
	            geom_point( aes(x=as.numeric(time), y= rescale(abruptness_abs_tw5_tp1, to=c(0, 1)), group=key), 
	                        shape=23, fill='blue', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(abruptness_abs_tw5_tp1, to=c(0, 1)), group=key), 
	                   color='blue', alpha=0.2, size=0.2 ) + 
	            facet_wrap( ~ replicate.id, scales='free') +
	            scale_fill_manual( values= color$asv[[1]] ) +
	            labs(x="Day", title=i, y= "Number of DNA coies",
	                 subtitle='red=relative, blue=absolute')+
	            guides(fill=guide_legend(title='ASV', nrow=2, reverse=TRUE))+ 
	            theme_minimal()+
	            theme(legend.position='') +
	            scale_y_continuous(sec.axis=sec_axis(~.* 1, name='Abruptness (time window 5, tp 1)')) 
    plot(ggarea)
    ## ========================================================== ##
	
}
dev.off()

saveRDS(all, 'Table/03_04_sample_info.rds')
head(all)
############################################################################

all <- readRDS('Table/03_04_sample_info.rds')

## ========================================================== ##
## -- Relationship of richness and abruptness

treat = unique(all[,c(5, 10)])

means <- c()
for(i in 1:nrow(treat)){#i=1
	
	## ====================================================== ##
	## -- Initianl state richness
	tmp <- Subset(all, treat, i)
	tmp2 <- tmp[as.numeric(tmp$time)%in%1:5, ]
	
	initialrichs <- c()
	for(l in colnames(taxa)[1:7][-2]){#l=colnames(taxa)[1]
		ts <- dlist[[as.character(treat[i,2])]][rownames(tmp2), ]
	
		a <- t(Taxa.mat(ts, taxa, l))
		initialrichs <- c(initialrichs, mean(apply(a, 1, function(x){ sum(x>0) } )))
		
	}
	names(initialrichs) <- colnames(taxa)[1:7][-2]
	
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	## -- Mean of richness
	meanrichs <- c()
	for(l in colnames(taxa)[1:7][-2]){#l=colnames(taxa)[1]
		ts <- dlist[[as.character(treat[i,2])]][rownames(tmp), ]
	
		a <- t(Taxa.mat(ts, taxa, l))
		meanrichs <- c(meanrichs, mean(apply(a, 1, function(x){ sum(x>0) } )))
		
	}
	
	names(meanrichs) <- colnames(taxa)[1:7][-2]
	
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	
	means <- rbind(means, c(unlist(treat[i,]), 
						   abruptness_rel_tw5_tp1 =max(tmp[,'abruptness_rel_tw5_tp1'], na.rm=TRUE ), 
						   abruptness_abs_tw5_tp1 =max(tmp[,'abruptness_abs_tw5_tp1'], na.rm=TRUE ),
						   initial= initialrichs, mean=meanrichs))
						   
	## ===================================================== ##				
}


## ========================================================= ##
df <- as.data.frame(means)
df[,-c(1,2)] <- apply(df[,-c(1,2)], 2, function(x){ as.numeric(as.vector(x))} )
lf <- gather(df, key, value, -c(1:4))
lf$key <- factor(lf$key, level=c(paste('initial', colnames(taxa)[1:7][-2], sep='.'), paste('mean', colnames(taxa)[1:7][-2], sep='.')))

lf$key <- gsub('initial.', 'Initail state richness of ', lf$key)
lf$key <- gsub('mean.', 'Mean richness of ', lf$key)

lf$key <- factor(lf$key, level=c("Initail state richness of ID", "Initail state richness of Phylum", "Initail state richness of Class", 
"Initail state richness of Order",  "Initail state richness of Family", "Initail state richness of Genus", 
"Mean richness of ID",  "Mean richness of Phylum", "Mean richness of Class",         
"Mean richness of Order", "Mean richness of Family","Mean richness of Genus"))

lf$treat1 <- factor(as.character(lf$treat1), levels=c("Soil/medium A", "Water/medium A", "Soil/medium B", 
					"Water/medium B", "Soil/medium C", "Water/medium C"))
					
g1 <- ggplot(lf)+
	  geom_point(aes(x=value, y= abruptness_abs_tw5_tp1, color=treat1))+
	  scale_color_manual(values=color$treat[levels(unique(lf$treat1))], 
	  					 guide=guide_legend(title='', ncol=3, direction='horizontal'))+
	  theme_minimal()+
	  theme(panel.border=element_rect(fill=NA, size=0.8),
	  		legend.position='bottom')+
	  facet_wrap(~key, scales='free', ncol=3)+
	  labs(y='Maximum abundance abruptness (time window 5, tp 1)')
	  
g2 <- ggplot(lf)+
	  geom_point(aes(x=value, y= abruptness_rel_tw5_tp1, color=treat1))+
	  scale_color_manual(values=color$treat[levels(unique(lf$treat1))], 
	  					 guide=guide_legend(title='', ncol=3, direction='horizontal'))+
	  theme_minimal()+
	  theme(panel.border=element_rect(fill=NA, size=0.8),
	  		legend.position='bottom')+
	  facet_wrap(~key, scales='free', ncol=3)+
	  labs(y='Maximum composition abruptness (time window 5, tp 1)')	
	   	  
pdf(sprintf('%s/maximum_abruptness_richness.pdf', dir$figdir), w=12, h=15)
plot(g1); plot(g2); dev.off()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
g1 <- ggplot(all)+
	  geom_point(aes(x=as.numeric(time), y= beta.q0_against_replicate, color= as.factor(replicate.id)))+
	  geom_line(aes(x=as.numeric(time), y= beta.q0_against_replicate, color= as.factor(replicate.id)))+
	  scale_color_manual(values=color$replicate, 
	  					 guide=guide_legend(title='Replicate', ncol=4, direction='horizontal', byrow=TRUE))+
	  theme_minimal()+
	  theme(panel.border=element_rect(fill=NA, size=0.8),
	  		legend.position='bottom')+
	  facet_wrap(~treat1, scales='free', ncol=3)+
	  labs(y='Mean of beta diversity against replicate', x='Time')

g2 <- ggplot(all)+
	  geom_point(aes(x=as.numeric(time), y= rel.cvs, color= as.factor(replicate.id)))+
	  geom_line(aes(x=as.numeric(time), y= rel.cvs, color= as.factor(replicate.id)))+
	  scale_color_manual(values=color$replicate, 
	  					 guide=guide_legend(title='Replicate', ncol=4, direction='horizontal', byrow=TRUE))+
	  theme_minimal()+
	  theme(panel.border=element_rect(fill=NA, size=0.8),
	  		legend.position='bottom')+
	  facet_wrap(~treat1, scales='free', ncol=3)+
	  labs(y='Mean CV of each ASV dynamics (relative)', x='Time')

g3 <- ggplot(all)+
	  geom_point(aes(x=as.numeric(time), y= ta.cv, color= as.factor(replicate.id)))+
	  geom_line(aes(x=as.numeric(time), y= ta.cv, color= as.factor(replicate.id)))+
	  scale_color_manual(values=color$replicate, 
	  					 guide=guide_legend(title='Replicate', ncol=4, direction='horizontal', byrow=TRUE))+
	  theme_minimal()+
	  theme(panel.border=element_rect(fill=NA, size=0.8),
	  		legend.position='bottom')+
	  facet_wrap(~treat1, scales='free', ncol=3)+
	  labs(y='CV of total abundace', x='Time')	

	  	    	  	 	  
pdf(sprintf('%s/Dynamics_indices.pdf', dir$figdir), w=15, h=9)
plot(g1); plot(g2); plot(g3); dev.off()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
summary(aov(abruptness_abs_tw5_tp1 ~treat1, data=all))
pairwise.t.test(all$abruptness_abs_tw5_tp1, all$treat1)	

g1 <- ggplot(all)+
	  geom_line(aes(x= time, y= abruptness_rel_tw5_tp1, group=replicate.id))+
	  facet_wrap(~treat1)
g2 <- ggplot(all)+
	  geom_line(aes(x= time, y= abruptness_abs_tw5_tp1, group=replicate.id))+
	  facet_wrap(~treat1)

	  
pdf(sprintf('%s/abruptness_size.pdf', dir$figdir), w=15,h=10)
plot(g1); plot(g2); dev.off()

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
plot(all$ewsPCA.PC1_tp10)
ggplot(all)+
	  geom_point(aes(x= ewsPCA.PC1_tp10, y= abruptness_rel_tw5_tp1, color= bray_against_replicate))+
	  facet_wrap(~treat1, scales='free')+scale_color_gradient2(midpoint=0.5)
		