############################################################################
####
#### R script for Fujita (2019)
####
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

#
# this script did not use for the paper (https://doi.org/10.1101/2022.08.23.505041)
#


## -- Loading Function and Library
source('functions/functions.R')
# -- Load library and functions
load.lib(c('vegan', 'e1071', 'doParallel'))

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/sample_info.rds')

time.window = c(10,15,20)

for.parallel(8)
#######################################################################################
## -- functions
normalize <- function(x){ (x-min(x)) / (max(x)-min(x)) }

cv <- function(x){ sd(x)/mean(x)}

#######################################################################################
ews <- c()

for(i in names(dlist)[-7]){ # i=names(dlist)[1]
    
    ## ========================================================== ##
  	## -- Extracting one treatment matrix
	
	smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]

    ## ========================================================== ##
    ## -- ASVs base Early warning signals
    ews.sum <- NULL
	for(tt in time.window){ # tt=10
		
		## ========================================================== ##
		## Absolute value

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
		
		asvews.abs <- cbind(cv=apply(cvlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}), 
						   skewness=apply(skewlist, 1, function(x){ mean(x[x>0], na.rm=TRUE)}) )
		rownames(asvews.abs) <- rownames(smlsub)
	    
		ews.tmp <- as.data.frame(asvews.abs[rownames(smlsub),])
		colnames(ews.tmp) <- paste(colnames(ews.tmp), '_tp', tt, sep='')
		
		ews.sum  <- cbind(ews.sum, as.matrix(ews.tmp))
	}	
    ## ========================================================== ##
    
    ews <- rbind(ews, ews.sum)
}

saveRDS(ews, 'Table/03_Early_warning_signals.rds')

#######################################################################################

