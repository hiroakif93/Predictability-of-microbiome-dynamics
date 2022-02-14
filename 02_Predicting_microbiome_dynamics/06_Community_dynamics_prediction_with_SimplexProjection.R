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

for.parallel(40)
############################################################################
## -- Prediction by using simplex projection

predlist <- statall <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ]

    ## ============================================== ##
    simpres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=1
        
        asv <- ts[,k]
	      absasv <- ts[,k]
        stat <- model <- c()
        for(r in 1:8){#r=2
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            E <- bestE_bidirect(asv, E_range = E.range, lib = libmat[-r,], criteria = "rmse", show_fig = F)
            simplex.out <- simplex(asv, E = E,  tp = tp, tau = 1,
                                   lib = libmat[-r,], pred=libmat[r,], 
                                   stats_only=FALSE)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rhormse <- sapply(simplex.out$model_output, function(x){ #x=simplex.out$model_output[[1]]
                
                n= sum(absasv[libmat[r,1]:libmat[r,2]]>0)
                y <- na.omit(x[,2:3])
                rmse <-  sqrt(sum((y[,1]-y [,2])^2)/length(y)) / (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                
                if(!all(x[,2]==0, na.rm=TRUE)){
                	cor=cor(x[,2], x[,3], use='pairwise.complete.obs',method=c('pearson'))
                }else{
                	cor=NA
                }
                c(rho= cor, rmse= rmse) })
            
            stat <- rbind(stat, cbind(asv=colnames(ts)[k], replicate=r, simplex.out[,c(1,3)], t(rhormse)))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_simplex_tp%s', colnames(ts)[k], tp)
            
            tmpmat <- matrix(NA, ncol=10, nrow=110, dimnames=list(libmat[r,1]:libmat[r,2], colname))
            for(l in 1:nrow(simplex.out)){ #l=1
                
                tmp <- simplex.out$model_output[[l]]
                rownames(tmp)[which(!is.na(tmp$time))] <- tmp$time[which(!is.na(tmp$time))]
                tmpmat[,l] <- tmp[rownames(tmpmat),'pred']
            }
            rownames(tmpmat) <- rowname
            model <- rbind(model, tmpmat)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        }
        
        list(stat, model)
    }
    
    ## ============================== ##
    
    statics <- do.call(rbind, simpres[,1])
    modelsimp <- c()
    for(j in tp){#j=1
    	
      modeltmp <- sapply(simpres[,2], function(x){  x[,j] })
      modeltmp.names <- sapply(simpres[,2], function(x){ colnames(x)[j] })
      colnames(modeltmp) <- modeltmp.names

    	modelsimp <- rbind(modelsimp, cbind(tp=j, modeltmp))
    }
	
    ## ============================== ##
    
    predlist[[i]] <- modelsimp
    statall[[i]] <- statics
    ## ============================== ##
    
}

saveRDS(predlist, 'Table/02_06_simplexModel.rds')    
saveRDS(statall, 'Table/02_06_simplexStatics.rds')    

    
