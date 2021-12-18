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
source("functions/BidirectSimplex_v0.7.4.R")
load.lib( c('rEDM'))

# -- Create directory to save
dir <- make.dir('06_community_prediction/06_05_output')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/03_matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/06_04_sample_info.rds')
color <- readRDS('Table/03_color_palette.rds')

# -- Decided regularized parameter range
time.window=c(1:6)
Erange=1:10

# -- Library point for SSR
libmat <- dlist[[7]]
lib <- cbind(dlist[[7]][,1], dlist[[7]][,1]+19)
pred.point <- dlist[[7]]

for.parallel(8)

############################################################################
## -- Function for simplex projection

simplex.f <- function(ts.data, Erange=1:8, libsize=20, tp=1){#ts.data=asv[1:110]; E=10; lib=libmat[-1,]; pred=libmat[1,]
    
    
    ## ================================= ##
    ## -- Preparing
   
    ## ================================= ##
    # Make an output dataframe
    output <- data.frame(index = 1:length(ts.data), obs=ts.data, pred = rep(NA, length(ts.data)) )
    
    s <- proc.time()[3]
    for(pp in libsize:(length(ts.data)-tp )+1 ){#pp =78
        
        if(!any(is.na(ts.data[pp]))){
        	
        	ts.lib = ts.data[(pp-libsize+tp):pp-tp]
        	Ebest=order(simplex(ts.data[(pp-libsize+tp):pp-tp], E= Erange, stats_only=FALSE)$rmse)[1]
        	
        	# -- Number of nearest neighbor points
    		num_neighbors= Ebest +1

    		# -- Time-delayed matrix
    		ts.embed = embed(c( rep(NA, Ebest-1), ts.lib), Ebest)

            set.target <- matrix( ts.embed[nrow(ts.embed)-(tp-1),], nrow=nrow(ts.embed), ncol=ncol(ts.embed), byrow=TRUE )
            dist <- sqrt( rowSums( (ts.embed - set.target)^2 ))
            
            nn <- order(dist)[-1][1:num_neighbors]
            
            min.distance <- dist[nn[1]]
            
            if(min.distance == 0){
                weights <- rep.int(0.000001, times = num_neighbors)
                weights[dist[nn] == 0] <- 1
            }else{
                weights <- exp(-dist[nn]/min.distance)
                weights[weights < 0.000001] <- 0.000001
            }
            total.weight <- sum(weights)
            output[output$index==(pp), 'pred'] <- (weights %*% ts.lib[(nn+tp)]) / total.weight
        }
        
    }
    
    ## ================================= ##
    return(output)
}

############################################################################
## -- Prediction by using simplex projection

error.sum <- predict.comm <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[2]

    ## ========================================================== ##
    ## -- Extracting one treatment matrix
	ts <- dlist[[i]]	
	smlsub <- sml[sml$treat1==i, ]
	
	## ========================================================== ##
	predmat <- foreach(j=1:ncol(ts))%dopar%{ #j=2
		
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
	    asv <- ts[,j]
		predmerge <- c()
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
		for(r in 1:8){#r=1
			
			mat <- matrix(NA, nrow=110, ncol=max(time.window))
			for(t in 0:(109-lib[1,2])){#t=0
				
			    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
				## -- Optimal E for each time ahead forecasting
			    asv2=asv[lib[r, 1]:lib[r, 2]+t]
				Etmp <- simplex(asv2, lib=c(1,19), 
								E=Erange, tau=1, tp= time.window, stats_only=F)
				
				## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
				reppred <- sapply(time.window, function(tt){ #tt=5
					Ebset <- order(Etmp[ Etmp$tp==tt, 'rmse'])[1]
					simp <- simplex(asv2, lib=c(1,19), pred=c(1, 20),
									E= Ebset, tau=1, tp= tt, stats_only=F)
					model <- simp$model_output[[1]]	
					model[which(model$time== 20), 'pred']
					
				})
				
				mat[lib[1, 2]+t+1, ] <- reppred
			}
			predmerge <- rbind(predmerge, mat)
		}
		
		return(predmerge)
	}
	
	## ======================================== ##
	
	predall <- c()
	for(tt in time.window){
		pred <- sapply(predmat, function(x) x[,tt])
		predall <- rbind(predall, cbind(tp=tt, pred))
	}
	
	## ======================================== ##

	predict.comm[[i]] <- predall
	
}

saveRDS(predict.comm, 'Table/06_06_prediction_comunity.rds')

############################################################################
# all <- cbind(infolist, error.sum[rownames(infolist), ])
predict.comm <- readRDS('Table/06_06_prediction_comunity.rds')

############################################################################
statsumm <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[2]
	
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
    
	ts <- dlist[[i]]
	predmat <- predict.comm[[i]]
	smlsub <- sml[sml$treat1==i, ]
	
	## ========================================================== ##
	
	direction <- euclid <- euclid2 <- obsdis <- preddis <- c()
	for(tt in time.window){#tt=5
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    ## -- No lag matrix
	    obs <- na.omit(ts); obs <- obs/rowSums(obs)
		pred <- predmat[predmat[,1]==tt, -1]; pred <- pred/rowSums(pred)
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Lag matrix
		dlyobs <- dlypred <- obs
		dlyobs[] <- dlypred[] <- NA
		
		for(r in 1:8){ #r=2
		    
		    obssub <- obs[(1:110)+110*(r-1), ]
		    prdsub <- pred[(1:110)+110*(r-1), ]
		    
		    dlyobs[((tt+1):110)+110*(r-1), ] <- obssub[1:(110-tt),]
		    dlypred[((tt+1):110)+110*(r-1), ] <- prdsub[1:(110-tt),]
		}

		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- (Xt - Xt-tau) - (Yt - Xt_tau)
		any(is.na(dlyobs[6,]))
		
		obs.bray <- rowSums(abs(obs-dlyobs)) / ( rowSums(obs) + rowSums(dlyobs) )
		pred.bray <- rowSums(abs(pred-dlyobs)) / ( rowSums(pred) + rowSums(dlyobs) )
		obs.eu <- sqrt( rowSums( (obs - dlyobs)^2 ))
		pre.eu <- sqrt( rowSums( (pred - dlyobs)^2 ))
		
		direction <- cbind(direction, abs(obs.bray-pred.bray))
		euclid <- cbind(euclid, abs(obs.eu-pre.eu))
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- (Xt - Xt-tau) - (Yt - Yt_tau)
		obs.eu <- sqrt( rowSums( (obs - dlyobs)^2 )) 
		pre.eu <- sqrt( rowSums( (pred - dlypred)^2 ))
		
		euclid2 <- cbind(euclid2, obs.eu-pre.eu)
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
		obsdis <- cbind(obsdis, obs.eu)
		preddis <- cbind(preddis,pre.eu)
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
	}
	## ========================================================== ##
	colnames(euclid) <- paste('euclid_t-tau_', 'tp', time.window, sep='')
	colnames(euclid2) <- paste('distance_obs-pred_', 'tp', time.window, sep='')
	colnames(direction) <- paste('braycurtis_t-tau_', 'tp', time.window, sep='')
	colnames(obsdis) <- paste('obs_euclid_', 'tp', time.window, sep='')
	colnames(preddis) <- paste('pred_euclid_', 'tp', time.window, sep='')
	
	statsumm <- rbind(statsumm, cbind(smlsub,  euclid, euclid2, obsdis, preddis))	
}

saveRDS(statsumm, 'Table/06_06_sample_info.rds')
