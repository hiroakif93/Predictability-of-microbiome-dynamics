############################################################################
####
#### R script for Fujita (2019)
####
#### CCM
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################
ran.seed <- 1234
set.seed(ran.seed)

## -- Loading Function and Library
source('functions/functions.R')
source("functions/BidirectSimplex_v0.7.4.R")
source('functions/Extended_Smap_v1.R')
source('functions/Extended_SSR_v1.R')

load.lib( c("rEDM", "glmnet", "tidyr", "doParallel") )

# -- Create directory to save
dir <- make.dir('03_Warning signals of critical transitions/Stablity')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat <- dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

# -- CCM threshold 
tw.th <- 3

time.horizon <- 1:4
Erange <- 1:10
srrogate.ittr <- 2500

# -- Decided regularized parameter range
alp = 0 # if 1 is lasso, 0 is ridge,
theta_test <-  c(0, 0.01,0.05,0.1,0.2,0.5,1,2,4,8)
lambda_test <- c(0, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2)
theta.lambda <- expand.grid(theta_test, lambda_test)

for.parallel(8)
############################################################################

statall  <- c()

for(i in names(dlist)[-7]){ # i=names(dlist)[2]
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
    ts <- dlist[[i]]
    smlsub <- sml[which(sml$treat1==i),]
    
    ## ========================================================== ##
	## -- Normalize
    ts.scale <- c()
    for(l in 1:8){#l=1
    	
    	tstmp <- ts[rownames(smlsub[smlsub$replicate.id==l,]),]
    	ts.scale <- rbind(ts.scale, rbind(as.matrix( apply(tstmp, 2, scale_2)), NA) )
    	
    }
    rownames(ts.scale) <- rownames(ts)
	
	############################################################################
	############################################################################
	## -- Data compile
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	## -- Combination list
	comb.tmp <- expand.grid(colnames(ts.scale), colnames(ts.scale))
	comb.tmp <- as.matrix(comb.tmp[which(comb.tmp[,1]!= comb.tmp[,2]), ])
	co.occur <- apply(comb.tmp, 1, function(x){ # x= comb.tmp[1,]
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- co-occure filtering
		y <- embed( c(rep(NA, tw.th-1), ts[, x[1]]), tw.th)
		z <- embed( c(rep(NA,tw.th-1), ts[, x[2]]), tw.th)
        
        appered.tw <- cbind(rowSums(y>0), rowSums(z>0))
        
        tmp <- c()
        for(r in 1:nrow(libmat)){#r=1
        	
        	size = c(libmat[r,1]:libmat[r,2])
        	tmp <- c(tmp, sum(rowSums( appered.tw[size, ]>2 )==2, na.rm=TRUE))
        }
        
        sum(tmp>0)
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    })
    
    ts.comb <- cbind(comb.tmp[ which(co.occur>2), ], co.occur=co.occur[ which(co.occur>2)]); dim(ts.comb)
    
    ############################################################################
    ############################################################################
    ## -- Causality Test by CCM
    ## test ts.comb <- ts.comb[sample(1:nrow(ts.comb), 24), ]
    
    ## ========================================================== ##
	replicateres <- c()
	for(r in 1:8){#r=1
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Detremining embedding dimension
        Elist <- foreach(e=colnames(ts.scale), .combine=rbind, .packages='rEDM') %dopar% { # e=colnames(ts)[2]
            
            simplex.res <- simplex(ts.scale[,e], lib=libmat[-r,], pred=libmat[-r,], E=Erange, silent = TRUE)
            order( simplex.res$rmse )[1]
        }
        Elist2 <- data.frame(Elist, row.names=colnames(ts.scale))
        
        ccm.comb <- cbind(ts.comb, E=Elist2[ts.comb[,1], ] )
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ccm.result <- foreach( rows=1:nrow(ccm.comb), .combine=rbind, .packages='rEDM') %dopar% { # rows=1	
            
            ## -- Targert variables
            ## lib.var is a variable used for library, and predicted.var is a predicted variable 
            lib.var = as.character(ccm.comb[rows, 1])
            predicted.var = as.character(ccm.comb[rows, 2])
 
            lib.E = as.numeric(ccm.comb[rows, 4] )
            
            ts.sub <- ts[libmat[r, 1]:libmat[r,2], ]
            cooccur <- sum( ts.sub[,lib.var]>0 & ts.sub[,predicted.var]>0, na.rm=TRUE )
         
            ## -- Run ccm
            ccm.tmp <- c()
            for(tp in time.horizon){ # tp=time.horizon[1]
                ccm.res <- ccm( block=ts.scale[,c(lib.var, predicted.var)], 
                                lib=libmat[-r,], pred=libmat[-r,],
                                lib_column=1, target_column=2, lib_size=c(10, 750),
                                E=lib.E, tau=1, tp= tp,
                                silent=TRUE)
                
                ccm.mean <- aggregate(ccm.res[,-2:0 + ncol(ccm.res)], by=list(ccm.res$lib_size), mean, na.rm=TRUE)			
                
                base=data.frame(replicate=r, lib.var=lib.var, predicted.var= predicted.var, cooccur_in_target= cooccur, E=lib.E, tp=tp)
                delta=ccm.mean[2, -1] -ccm.mean[1, -1]
                
                ccm.tmp <- rbind(ccm.tmp, cbind(base, delta=delta, Lmin=ccm.mean[1, -1], Lmax=ccm.mean[2, -1]))
            }
            
            ccm.tmp$best <- order(ccm.tmp$delta.rmse)
            return(ccm.tmp)
            
        }
		
		replicateres <- rbind(replicateres, ccm.result)
    }   
    
	############################################################################
	############################################################################
    ## -- Surrogate test
    
    ## ========================================================== ##
    ## -- Two type of surrogate deta
        
    replicateres$pvalue_fourier <- NA
    replicateres$pvalue_seasonal <- NA
         
    for(surr in which(replicateres$best==1 & replicateres$delta.rmse<0 & 
    					replicateres$cooccur_in_target > 2)){ 
        #surr=which(replicateres$best==1 & replicateres$delta.rmse<0 & replicateres$cooccur_in_target > 2)[1]
    	
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Extract parameter
        lib.var = as.character(replicateres[surr, 'lib.var'])
        predicted.var = as.character(replicateres[surr, 'predicted.var'])
        E = as.numeric(replicateres[surr, 'E'] )
        tp = as.numeric(replicateres[surr, 'tp'] )
		replicateID <- as.numeric(replicateres[surr, 'replicate'] )
        delta.rmse = as.numeric(replicateres[surr, 'delta.rmse'] )
        
		## -- Extract ASVs
		asv.lib <- ts.scale[,lib.var]
		asv.pred <- ts.scale[,predicted.var]
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
        ## -- Making surrogate data
        seasn.data <- forir.data <- matrix(NA, ncol= srrogate.ittr, nrow=nrow(ts.scale))
        for(r in 1:8){#r=1
        	
        	rows <- libmat[r,1]:libmat[r,2]
        	sub <- ts.scale[rows, predicted.var]
        	
        	seasn.data[rows,] <- make_surrogate_data( ts=sub, method='seasonal', num_surr = srrogate.ittr, T_period = 110)
        	forir.data[rows,] <- make_surrogate_data( ts=sub, method='ebisuzaki', num_surr = srrogate.ittr)
        	
        }
        
        ## -- Run CCM
        seasn.test.pvalue <- foreach(num=1:srrogate.ittr, .combine=c, .packages='rEDM') %dopar% { #num=1
                
                block <- cbind(asv.lib, seasn.data[,num])
                ccm.res <- ccm( block= block, 
                                lib= libmat[-replicateID,], pred= libmat[-replicateID,],
                                lib_column=1, target_column=2, lib_size=c(10, 100),
                                E=E, tau=1, tp= tp,
                                silent=TRUE)
                
                ccm.mean <- aggregate(ccm.res[,-2:0 + ncol(ccm.res)], by=list(ccm.res$lib_size), mean, na.rm=TRUE)		
                ccm.mean[2, 4] -ccm.mean[1, 4]		
        }
            
        replicateres[surr, 'pvalue_seasonal'] <- sum(seasn.test.pvalue < delta.rmse) / srrogate.ittr
        
        forir.test.pvalue <- foreach(num=1:srrogate.ittr, .combine=c, .packages='rEDM') %dopar% { #num=1
                
                block <- cbind(asv.lib, forir.data[,num])
                ccm.res <- ccm( block= block, 
                                lib= libmat[-replicateID,], pred= libmat[-replicateID,],
                                lib_column=1, target_column=2, lib_size=c(10, 100),
                                E=E, tau=1, tp= tp,
                                silent=TRUE)
                
                ccm.mean <- aggregate(ccm.res[,-2:0 + ncol(ccm.res)], by=list(ccm.res$lib_size), mean, na.rm=TRUE)		
                ccm.mean[2, 4] -ccm.mean[1, 4]		
        }
            
        replicateres[surr, 'pvalue_fourier'] <- sum(forir.test.pvalue < delta.rmse) / srrogate.ittr
            
    }
	############################################################################
	############################################################################
	
	result <- cbind(treatment=i, replicateres)       
	saveRDS(result, sprintf('%s/CCM_crosvaridation_%s.rds', dir$rdsdir, gsub('/', '_', i)))
	
	############################################################################
	############################################################################
	## -- Regularized multivarite S-map

	## ========================================================== ##
	ccmres <- result[which(result$pvalue_fourier<0.05), ]
	pattern <- unique(ccmres[, c(2,3)])
	
	coefficients<- foreach(ptn = 1:nrow(pattern), .combine=rbind, .packages=c("rEDM", "glmnet")) %dopar% { #ptn=50
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Target
        causeset <- Subset(ccmres, pattern, ptn)
        causeset <- causeset[order(causeset$delta.rmse),]
        
        lib.var <- as.character(unique(causeset[,'lib.var']))
        predicted.var <- as.character(unique(causeset[,'predicted.var']))
       
        target.rep = pattern[ptn,1]
        target.row = libmat[target.rep,1] : libmat[target.rep,2] 
        
        coef_tp <- stat_tp <- model_tp <- c()
        for(tp in 1:10){ #tp=1
        	
            ## ================================================= ##
            ## -- Make time series block
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##	
            ## --Cause
            causeblock <- causets <- as.matrix(ts.scale[, predicted.var])
            causeblock[] <- NA
            for(row in 1:nrow(causeset)){ # row=1
            	
                delay.tp <- rep(NA, (abs(causeset$tp[row])-1))
                
                for(n in 1:8){ #n=1
					librow <- libmat[n,1] : libmat[n,2]       	
					causeblock[librow, row] <- c(delay.tp, causets[librow,row])[1:110]
		        }
    
            }
            
            if( ncol(causeblock)>1){
                target <- which(apply(causeblock, 2, sd, na.rm=TRUE)>0)
                causeblock <- causeblock[,target]; cause <- names(target)
            }else{
                cause <- as.character(unique(causeset[,'predicted.var']))
            }
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## --Effect
            libts <- ts.scale[, lib.var]
            if( ncol(causeblock) >= unique(causeset[,'E']) ){
                libblock <- matrix(libts, ncol=1)
                lag <- 'lag_1'
            }else{
                dly <- unique(causeset[,'E']) - ncol(causeblock) 
                libblock <- embed( c(rep(NA, dly-1 ), libts), dly )
                
                lag <- paste('lag', 1:dly, sep='_')
            }
            
            ## ========================================================== ##
            # -- Decide optimal theta and lamda values by RMSE
            
            para <- c()
            for(cols in 1:ncol(libblock)){
                
                for(row in 1:nrow(theta.lambda)){
                    
                    block <- cbind(libblock[,1:cols], causeblock)
                    
                    tmp <- extended_lnlp(block, theta = theta.lambda[row,1], lambda = theta.lambda[row,2],
                                         lib = libmat[-target.rep, ],  pred=libmat[-target.rep, ],            
                                         method = "s-map", regularized = T, alpha = alp, 
                                         random_seed = ran.seed, no_parallel = T, tp=tp)
                    para <- rbind(para,data.frame('theta' = theta.lambda[row,1],'lambda' =  theta.lambda[row,2], 'cols'=cols, tmp $stats))
                    
                }
                
            }
            
            best.para <- para[para $rmse==min(para$rmse),]
            block <- cbind(libblock[,1:best.para$cols], causeblock)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## -- Infer S-map coefficient by using regularized multivariate S-map method
           	res.tmp <- extended_lnlp(block, lib=libmat[-target.rep, ],  pred=libmat[target.rep, ],
                                     method = "s-map", save_smap_coefficients = T, 
                                     theta = best.para$theta, lambda = best.para$lambda, regularized = T,
                                     random_seed = ran.seed, no_parallel = F,tp= tp,
                                     alpha = alp)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # -- Store results
           	# -- model
           	model <- data.frame(treat=i, effect=as.character(unique(causeset[,'lib.var'])),
           	                    cbind(tp=tp, replicate=target.rep, time=1:110,res.tmp[[1]][target.row,-1]))
           	model_tp <- rbind(model_tp, model)
           	# -- statics
           	stat <- data.frame(treat=i,effect=as.character(unique(causeset[,'lib.var'])), tp=tp, replicate=target.rep,  as.matrix(res.tmp[[2]]))
           	stat_tp <- rbind(stat_tp, stat)
           	
            # -- Coefficient
            coef <- as.matrix(res.tmp$smap_coefficients)[target.row,]
            
            coef.t <- as.data.frame( t(coef) )
            coef.t2 <- cbind(matrix(NA, ncol=10, nrow=nrow(coef.t)), coef.t[,-c(1:10, ncol(coef.t))], NA)
            
            rownames(coef.t2)[nrow(coef.t2)] <- 'intercept'
            rownames(coef.t2)[grep('c_',rownames(coef.t2))] <- c(lag[1:best.para$cols], cause)      
    
            tmp <- cbind( treat=i, tp=tp, replicate=target.rep, effect=as.character(unique(causeset[,'lib.var'])), cause=rownames(coef.t2)[-1], coef.t2[-1,])
            
            coef_tp <- rbind(coef_tp, tmp)
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            
        }
        
        return(list(coef_tp, stat_tp, model_tp))
    }
	
	saveRDS(coefficients, sprintf('%s/RMsmap_coef_v1_%s_list.rds',dir$rdsdir,  gsub('/','_',i)))
   
	############################################################################
	############################################################################
    statall[[i]] <- coefficients
	
}    
saveRDS(statall, 'Table/03_01_RMsmap_coef_list.rds')
       
        