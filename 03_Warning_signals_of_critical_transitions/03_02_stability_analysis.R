############################################################################
####
#### R script for Fujita (2019)
####
#### Interaction network analysis
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS')
#### 
############################################################################

ran.seed <- 123
set.seed(ran.seed)

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c("ggplot2", "tidyr", "scales"))

# -- Create directory to save
dir <- make.dir('03_Warning_signals_of_critical_transitions/Stability')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat <- dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

for.parallel(8)
########################################################################

fl <- list.files(dir$rdsdir, full.names=TRUE)
rmsmap <- lapply(fl[grep('RMsmap_coef_v1', fl)], readRDS)
is0710 <- c()
for(i in 1:length(rmsmap)) {tmp <- do.call(rbind, rmsmap[[i]][,1] ) ; is0710 <- rbind(is0710, tmp)}
########################################################################
treat <- (unique(is0710[,c(1,3)]))

statsum <- c()
eigenmat <- c()
for(i in 1:nrow(treat)){ #i=1
	
	## ========================================================== ##
	## -- Extracting one treatment matrix
		
	smlsub <- sml[sml$treat1==as.character(treat[i,1]) & sml$replicate.id ==as.character(treat[i,2]), ]
	ts <- dlist[[ as.character(treat[i,1]) ]][rownames(smlsub), ]
	    
	smltmp <- smlsub
	
	for(tp in 1){#tp=1
		## ========================================================== ##
	    ## -- Extracting one treatment matrix
		
	    is <- Subset(is0710, treat, i)
		is <- is[is$tp==tp, ]	
		
	    ## ========================================================== ##
	    ## -- Deleating 'intercept' rows
	    issub <- is[ -grep('intercept', is[,'cause']), ]
	
		## ========================================================== ##
		issub2 <- issub
		issub2[,4:5] <- apply(issub2[,4:5], 2, as.character)
		issub2[grep('lag_1', issub2[,5]), 5] <- issub2[grep('lag_1', issub2[,5]), 4]
		issub2[grep('lag', issub2[,5]), 5] <- apply(issub2[grep('lag', issub2[,5]), 4:5], 1, paste, collapse='_')
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		statmat <- matrix(NA, ncol=2, nrow=110, dimnames=list(NULL, paste(c('DS_0710', 'LSS_0710'), '_tp', tp, sep='')))
		statmat2 <- matrix(NA, ncol=8, nrow=110, 
		dimnames=list(NULL,  paste('is_0710_', c('mean', 'mean_abs', 'max', 'min', 'positive_is', 'negative_is', 'positive', 'negative'),'_tp', tp, sep='') ))
		
		causes <- unique(issub2$cause)
		causes <- causes[-grep('lag', causes)]
		eigenmat.tmp <- matrix(0, nrow=110, ncol=length(causes), dimnames=list(NULL, causes))
		
		for(p in 6:ncol(issub)){ #p=60
			
			## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
			## -- Interaction each day
			present <- colnames( ts )[which(ts[p-5,]>0)]
			tmp1 <- issub2[,c(1:5, p)]
			tmp2 <- tmp1[which(tmp1[,4] %in% present),]
			istmp <- tmp2[c(which(tmp2[,5] %in% present), grep('lag', tmp2[,5])),]
			
			## ============================================= ##
		    ## -- Interaction strength from X
		    isX <- istmp[-grep('lag_', istmp$cause), ]	
			ismat <- t(isX[,-c(1:5)])
			
			stat <- t(apply(ismat, 1, function(x){ #x=ismat[25, ]
						c(mean(x), mean(abs(x)), max(x), min(x), sum(x>0)/length(x), sum(x<0)/length(x), mean(x[x>0]), mean(x[x<0]) ) }))
			statmat2[p-5, ] <- stat

			## ============================================= ##
		
			## -- Making template Jacobian matrix
		    rowcol <- unique( c(as.character(istmp[,4]), as.character(istmp[,5])) )
		    
		    lag <- max(as.numeric(gsub('.*_lag_', '',rowcol[grep('lag', rowcol)])))
		    if(lag>0){
		    	asvname <- rowcol[-grep('lag', rowcol)]
		    }else{
		    	asvname <- rowcol
		    }
		    
		    mattmp <- matrix(0, ncol=length(asvname), nrow=length(asvname), 
		    				dimnames=list(asvname, asvname))
			
			## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
			# -- Binding time-lagged Jacobian matrix
			lags <- ifelse(lag-1>0, 1:(lag-1), 0)
			if(lag-1>0){
				for(l in 1:(lag-1)){
					
					lagname <- paste(asvname, sprintf('_lag_%s', (l+1)), sep='')
					tmp <- matrix(0, ncol=length(asvname), nrow=length(asvname), 
			    				dimnames=list(lagname, lagname))
			    				
					mattmp <- cbind(mattmp, tmp)
				}
				
				I <- diag(nrow(mattmp)* (lag-1) )
				O   <- matrix(0, nrow = nrow(mattmp) * (lag-1), ncol = nrow(mattmp))
				
				matA <- rbind(mattmp, cbind(I, O))
			}else{
				matA <- mattmp
			}	
			## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
			
			if( any( !is.na(issub[,p]) ) & nrow(matA)>0 ){
				
				## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
				## -- Insert interaction strength to Jacobian matrix
				
				for(n in 1:nrow(istmp)){ #n=33
					matA[istmp[n,4], istmp[n,5]] <- istmp[n,6]
				}
				
				## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##	
				## -- eigen value and trace of jacobian
				ev <- eigen(matA)
				maxev.p <- which.max( abs(ev$values) )
				max.ev <- max(abs(ev$values))
				
				evvec <- ev$vectors[,maxev.p]
				evvec.df <- data.frame(asv=sapply(strsplit(colnames(matA), '_lag'), '[', 1), evvec)	
				cont <- aggregate(evvec.df[,2], by=list(evvec.df$asv),sum); rownames(cont) <- cont[,1]
				
				eigenmat.tmp[p-5, ] <- Re(cont[colnames(eigenmat.tmp),2])
				
				#re.ev  <- Re(ev)
				#Max_re <- max(abs(re.ev))
				#max.ev <- abs(Re(ev[abs(Re(ev)) == Max_re]))
				if (length(max.ev) == 1) eigen_tmp <- max.ev
				if (length(max.ev) != 1) eigen_tmp <- max.ev[1]
			
				trJ <- sum(diag(matA))
				
				statmat[p-5, ] <- c(eigen_tmp, trJ)
			}
			
		}
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		smltmp <- cbind(smltmp, statmat, statmat2)
	}
    
	statsum <- rbind(statsum, smltmp)
	
	rownames(eigenmat.tmp) <- rownames(smltmp)
	eigenmat[[i]] <- eigenmat.tmp
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
}

saveRDS(statsum, 'Table/03_Stability.rds')
saveRDS(list(treat=treat, eigenmat), 'Table/03_eigenVector.rds')

## ========================================================== ##

