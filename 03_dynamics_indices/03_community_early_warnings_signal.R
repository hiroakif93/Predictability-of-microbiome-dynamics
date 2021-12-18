############################################################################
####
#### R script for Fujita (2019)
####
#### Unifraq
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################


## -- Loading Function and Library
source('functions/functions.R')
# -- Load library and functions
load.lib(c('scales', 'ggplot2', 'RColorBrewer', 'cowplot', 'vegan', 'e1071', 'tidyr', 'doParallel'))

# -- Create directory to save
dir <- make.dir('03_dynamics_indices/03_EWS')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist <- readRDS("Table/03_matrixList.rds")
sml  <- readRDS("Table/smaple_info_03_02.rds")[,-c(14,16,17)]
color <- readRDS('Table/03_color_palette.rds')

time.window = c(10,15,20)

for.parallel(8)
#######################################################################################
## -- functions
normalize <- function(x){ (x-min(x)) / (max(x)-min(x)) }

cv <- function(x){ sd(x)/mean(x)}

#######################################################################################
ews <- infoall <- c()

for(i in names(dlist)[-7]){ # i=names(dlist)[1]
    
    ## ========================================================== ##
  	## -- Extracting one treatment matrix
	
	smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]
	rel.ts <- ts/rowSums(ts)
		 
    ## ========================================================== ##
    ## -- ASVs base Early warning signals
    ews.sum <- NULL
	for(tt in time.window){ # tt=10
		
		## ========================================================== ##
		## Absolute value
		arlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=2
			
			ews.rep <- lapply(1:8, function(r){#r=1
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- scale_2(ts[rownames(sample.rep), cols])

				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				ar <- apply(dly.asv, 1, function(x){ ifelse( any(is.na(x)) | sd(x, na.rm=TRUE) == 0, 0, ar(x, order=1)$ar) })

				return(ar)
			})
			
			return(do.call(c, ews.rep))
		 }
		
		cvlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=1
			
			ews.rep <- lapply(1:8, function(r){
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- ts[rownames(sample.rep), cols]
				
				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				CV <- apply(dly.asv, 1, cv)
				
				return(CV)
			})
			
			return(do.call(c, ews.rep))
		}
		
		skewlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=1
			
			ews.rep <- lapply(1:8, function(r){
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- ts[rownames(sample.rep), cols]
				
				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				skew <- apply(dly.asv, 1, skewness)

				return(skew)
			})
			
			return(do.call(c, ews.rep))
		}
		
		covlist <- matrix(NA, ncol=1, nrow=880, dimnames=list(rownames(smlsub), NULL))
		for(r in 1:8){ #r=1
			sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
			tssub <- ts[rownames(sample.rep),]
			tssub <- apply(tssub, 2, scale_2)
			for(tp in tt:110){#tp=10
				covlist[rownames(sample.rep)[tp],] <- max(eigen(var(tssub[(tp-(tt-1)):tp,]))$value)
			}
		
		}
		
		asvews.abs <- cbind( ar1=apply(arlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}), 
						 cv=apply(cvlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}), 
						 skewness=apply(skewlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}), cov= covlist)
		rownames(asvews.abs) <- rownames(smlsub)

		## ========================================================== ##
		## Relative value
		arlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=1
			
			ews.rep <- lapply(1:8, function(r){#r=1
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- rel.ts[rownames(sample.rep), cols]
				
				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				ar <- apply(dly.asv, 1, function(x){ ifelse( any(is.na(x)) | sd(x, na.rm=TRUE) == 0, 0, ar(x, order=1)$ar) })

				return(ar)
			})
			
			return(do.call(c, ews.rep))
		 }
		
		cvlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=1
			
			ews.rep <- lapply(1:8, function(r){
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- rel.ts[rownames(sample.rep), cols]
				
				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				CV <- apply(dly.asv, 1, cv)
				
				return(CV)
			})
			
			return(do.call(c, ews.rep))
		}
		
		skewlist <- foreach(cols=1:ncol(ts), .combine=cbind)%dopar%{ #cols=1
			
			ews.rep <- lapply(1:8, function(r){
				sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
				asv <- rel.ts[rownames(sample.rep), cols]
				
				dly.asv <- embed( c(rep(NA, tt-1), asv), tt)
				skew <- apply(dly.asv, 1, skewness)

				return(skew)
			})
			
			return(do.call(c, ews.rep))
		}
		
		covlist <- matrix(NA, ncol=1, nrow=880, dimnames=list(rownames(smlsub), NULL))
		for(r in 1:8){ #r=1
			sample.rep <- smlsub[ which(smlsub$replicate.id==r), ]
			relsub <- rel.ts[rownames(sample.rep),]
			
			for(tp in tt:110){#tp=10
				covlist[rownames(sample.rep)[tp],] <- max(eigen(var(relsub[(tp-(tt-1)):tp,]))$value)
			}
		
		}
		
		
		asvews.rel <- cbind(ar1=apply(arlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}), 
							cv=apply(cvlist, 1, mean, na.rm=TRUE), skewness=apply(skewlist, 1, mean, na.rm=TRUE), cov= covlist)
		rownames(asvews.rel) <- rownames(smlsub)

		## ========================================================== ##
    	## -- diversity or total abundance base Early warning signals
		a.div <- apply( ts, 1, function(x){ sum(x>0) })
		a.div.dly <- embed(c( rep(NA, tt-1), a.div), tt)
		div.ews <- t(apply(a.div.dly, 1, function(x){
			c( ar1= ifelse(any(is.na(x)) | sd(x, na.rm=TRUE)==0, 0, ar(x, order=1)$ar), cv=cv(x), skewness=skewness(x))
		}))

		total <- apply( ts, 1, sum)
		total.dly <- embed(c( rep(NA, tt-1), total), tt)
		total.ews <- t(apply(total.dly, 1, function(x){
			c( ar1= ifelse(any(is.na(x)) | sd(x, na.rm=TRUE)==0, 0, ar(x, order=1)$ar), cv=cv(x), skewness=skewness(x))
		}))
		rownames(div.ews) <- rownames(total.ews) <- rownames(ts)
		colnames(div.ews) <- paste('alpha_',colnames(div.ews), sep='' )
		colnames(total.ews) <- paste('total.abun_',colnames(total.ews), sep='')

		ews.tmp <- cbind(asv.abs=as.data.frame(asvews.abs[rownames(smlsub),]),asv.rel=asvews.rel[rownames(smlsub),], 
						 #ewsPCA.abs=as.data.frame(abs.pca$x)[rownames(smlsub), 1:2],
						 #ewsPCA.rel=as.data.frame(rel.pca $x)[rownames(smlsub), 1:2],
						 div.ews[rownames(smlsub),], total.ews[rownames(smlsub),])
		pca <- as.data.frame(prcomp( na.omit(ews.tmp[,-grep('ar', colnames(ews.tmp))]))$x)[rownames(smlsub),1:2]				 
		ews.tmp <- cbind(ews.tmp, ewsPCA=pca)
		colnames(ews.tmp) <- paste(colnames(ews.tmp), '_tp', tt, sep='')
		
	
		ews.sum  <- cbind(ews.sum, as.matrix(ews.tmp))
	}	

    ews <- rbind(ews, ews.sum)
    ## ========================================================== ##

    smlsub2 <- cbind(smlsub, ews.sum[rownames(smlsub),])
    lfarea <- gather( cbind(smlsub2, #tmp.domi,
                            ts[rownames(smlsub2),]), key, value, -c(1:(ncol(smlsub2)))) 
    lfarea$key <- factor( lfarea$key, levels= names(sort(colSums(ts, na.rm=TRUE))) )
    ggarea1 <- ggplot(lfarea)+
	            geom_area( aes(x=as.numeric(time), y=value, group=key, fill=key ))+
	            geom_point( aes(x=as.numeric(time), y= rescale(ewsPCA.PC1_tp10, to=c(0, max(ts, na.rm=TRUE))), group=key), 
	                        shape=23, fill='red', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(ewsPCA.PC1_tp10, to=c(0, max(ts, na.rm=TRUE))), group=key), 
	                   color='red', alpha=0.2, size=0.2 ) + 
                geom_point( aes(x=as.numeric(time), y= rescale(ewsPCA.PC2_tp10, to=c(0, max(ts, na.rm=TRUE))), group=key), 
                            shape=23, fill='blue', size=1 ) +	
                geom_line( aes(x=as.numeric(time), y= rescale(ewsPCA.PC2_tp10, to=c(0, max(ts, na.rm=TRUE))), group=key), 
                           color='blue', alpha=0.2, size=0.2 ) +
	            facet_wrap( ~ replicate.id, scales='free') +
	            scale_fill_manual( values= color$asv$ID[colnames(ts)] ) +
	            labs(x="Day", title=i, y= "Number of DNA coies",
	                 subtitle='red=asv base cv, blue=asv base skewness')+
	            guides(fill=guide_legend(title='ASV', nrow=2, reverse=TRUE))+ 
	            theme_minimal()+
	            theme(legend.position='')   
     ggarea2 <- ggplot(lfarea)+
	            geom_area( aes(x=as.numeric(time), y=value, group=key, fill=key ), position='fill')+
	            geom_point( aes(x=as.numeric(time), y= rescale(ewsPCA.PC1_tp10, to=c(0, 1)), group=key), 
	                        shape=23, fill='red', size=1 ) +	
	            geom_line( aes(x=as.numeric(time), y= rescale(ewsPCA.PC1_tp10, to=c(0, 1)), group=key), 
	                   color='red', alpha=0.2, size=0.2 ) + 
                geom_point( aes(x=as.numeric(time), y= rescale(ewsPCA.PC2_tp10, to=c(0, 1)), group=key), 
                            shape=23, fill='blue', size=1 ) +	
                geom_line( aes(x=as.numeric(time), y= rescale(ewsPCA.PC2_tp10, to=c(0, 1)), group=key), 
                           color='blue', alpha=0.2, size=0.2 ) +
	            facet_wrap( ~ replicate.id, scales='free') +
	            scale_fill_manual( values= color$asv$ID[colnames(ts)] ) +
	            labs(x="Day", title=i, y= "Number of DNA coies",
	                 subtitle='red=asv base cv, blue=asv base skewness')+
	            guides(fill=guide_legend(title='ASV', nrow=2, reverse=TRUE))+ 
	            theme_minimal()+
	            theme(legend.position='')
          
    pdf(sprintf('%s/timeseries_abrupt_%s_abs.pdf', dir$figdir, gsub('/', '-', i)), w=15,h=10)
    plot(ggarea1);    dev.off()
    pdf(sprintf('%s/timeseries_abrupt_%s_rel.pdf', dir$figdir, gsub('/', '-', i)), w=15,h=10)
    plot(ggarea2);    dev.off()

    ## ========================================================== ##
    
    infoall <- rbind(infoall, smlsub2)
}

saveRDS(ews, 'Table/03_03_ews.rds')
saveRDS(infoall, 'Table/03_03_sample_info.rds')

#######################################################################################

