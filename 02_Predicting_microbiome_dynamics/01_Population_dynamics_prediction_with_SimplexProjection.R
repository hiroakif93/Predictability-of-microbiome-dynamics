############################################################################
####
#### R script for Fujita (2019)
####
#### Prediction by using simplex projection
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
source("functions/BidirectSimplex_v0.7.4.R")
load.lib( c('rEDM'))

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

# -- Decided regularized parameter range
E.range <- 1:20
time.horizon <- 1:10
tp=c(1:10)

for.parallel(detectCores())
############################################################################

predlist <- statall <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[6]
    
    ## ============================================== ##
    ts <- dlist[[i]]
    
    smlsub <- sml[sml$treat1==i, ]

    tmp <- lapply(1:8, function(x){
        y <- apply(ts[libmat[x,1]:libmat[x,2],], 2, scale_2)
        rbind(y, NA)
    })
    #tsscale <- do.call(rbind, tmp)
    tsscale <- apply(ts, 2, scale_2)

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    ## -- NA index
    index <- 1:nrow(ts); names(index) <- rownames(ts)
    naindex <- index[rownames(smlsub[smlsub$missing=='N', ])]
    
    ## ============================================= ##

    simpres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=38
        
        asv <- tsscale[,k]
		    absasv <- ts[,k]
        stat <- model <- c()
        for(r in 1:8){#r=4
            
            smlrep <- smlsub[smlsub$replicate.id==r,]
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            E <- bestE_bidirect(asv, E_range = E.range, lib = libmat[-r,], criteria = "rmse", show_fig = F)
            simplex.out <- simplex(asv, E = E,  tp = tp, tau = 1,
                                   lib = libmat[-r,], pred=libmat[r,], 
                                   stats_only=FALSE)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rhormse <- sapply(simplex.out$model_output, function(x){ #x=simplex.out$model_output[[10]]
                
                nan <- which(smlrep$missing=='N')[which(which(smlrep$missing=='N')%in%c((max(E.range)+1):110))]
                index=libmat[r, 1]:libmat[r,2]
                x <- x[x$time%in% index[nan], ]
              
                n= sum(absasv[libmat[r,1]:libmat[r,2]][nan]>0)
                y <- x[,2:3]
                
                rmse <-  sqrt(sum((y[,1]-y [,2])^2)/nrow(y)) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                rmse.weight <-  sqrt(sum((y[,1]-y [,2])^2)/n) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                rmse.weight <- ifelse(is.infinite(rmse.weight), NA, rmse.weight)
                
                if(!all(x[,2]==0, na.rm=TRUE)){
                	cor=cor(y[,1], y[,2], use='pairwise.complete.obs',method=c('spearman'))
                }else{
                	cor=NA
                }
                c(rho= cor, rmse= rmse, rmse.weight= rmse.weight) })
            
            stat <- rbind(stat, cbind(asv=colnames(tsscale)[k], replicate=r, simplex.out[,c(1,3)], t(rhormse)))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_simplex_tp%s', colnames(tsscale)[k], tp)
            
            tmpmat <- matrix(NA, ncol=10, nrow=110, dimnames=list(rowname, colname))
            for(l in 1:nrow(simplex.out)){ #l=1
                
                tmp <- simplex.out$model_output[[l]]
                tmpmat[,l] <- tmp[,'pred']
            }
            model <- rbind(model, tmpmat)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        }
        
        list(stat, model)
    }
    
    ## ============================== ##
    
    statics <- do.call(rbind, simpres[,1])
    modelsimp <- c()
    for(j in tp){#j=1
    	
    	modeltmp <- sapply(simpres[,2], function(x){ #x = simpres[,2][[1]]
    		x[,j]
    	})
    	modelname <- sapply(simpres[,2], function(x){ #x = simpres[,2][[1]]
    	  colnames(x)[j]
    	})
    	colnames(modeltmp) <- modelname
    	modelsimp <- rbind(modelsimp, cbind(tp=j, modeltmp))
    }
	
    ## ============================== ##
    
    predlist[[i]] <- modelsimp
    statall[[i]] <- statics
    ## ============================== ##
    
}

saveRDS(predlist, 'Table/02_01_simplexModel.rds')    
saveRDS(statall, 'Table/02_01_simplexStatics.rds')    
    
    
