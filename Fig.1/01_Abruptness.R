############################################################################
####
#### R script for Fujita (2019)
####
#### Abruptness
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'vegetarian', 'scales'))

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/sample_info.rds')

time.step = c(1:10)
previous=5
timewindow=c(3,5,7)

for.parallel(8)
############################################################################

abruptness <- all <- c()

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
	
}

saveRDS(all, 'Table/Abruptness.rds')

############################################################################
