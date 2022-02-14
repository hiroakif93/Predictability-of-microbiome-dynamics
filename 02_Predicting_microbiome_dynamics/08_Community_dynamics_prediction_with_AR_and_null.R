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
load.lib( c('rEDM'))

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

# -- Decided regularized parameter range
tp=c(1:10)
tw=10

for.parallel(8)
############################################################################

predlist <- statall <- meanmodels <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ============================== ##
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ]

    ## ============================== ##
    arres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=19
        
        asv <- ts[,k]
		absasv <- ts[,k]
        stat <- model <- c()
        
        for(r in 1:8){#r=1
          
            ## ============================== ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_AR_tp%s', colnames(ts)[k], tp)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ar.pred <-  matrix(NA, ncol=length(tp), nrow=110, dimnames=list(rowname, colname)) ; 
			arts <- asv[c(libmat[r,1]: libmat[r,2])]
            
			for(t in (tw):110){ #t=17
				
				if(sd(arts[1:t])>0){
					ar.coef <- ar(arts[1:t], aic=TRUE)$ar
					yt=arts[t:1][1:length(ar.coef)]
					
	            	intercept <- mean((1-ar.coef)*mean(yt))
	                
	                for(time.ahead in tp){#time.ahead=1
	            		
	            		pred <- sum(ar.coef*(yt))+ intercept
	            		yt <- c(pred, yt)[1:length(ar.coef)]
	            		
	            		if( (t+time.ahead) < 110) ar.pred[t+time.ahead, time.ahead] <- pred
	                }
                }else{
                	ar.pred[t,] <- 0
                }
			}         		
            
            ## ============================== ##
          
            rhormse <- apply(ar.pred, 2, function(x){ #x=ar.pred[,1]
            	
                n= sum(absasv[libmat[r,1]:libmat[r,2]]>0)
                y <- na.omit(cbind(asv[libmat[r,1]:libmat[r,2]], x))
                rmse <-  sqrt(sum((y[,1]-y [,2])^2)/length(y)) / (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                
                if(!all(y[,1]==0, na.rm=TRUE)){
                	cor=cor(y[,1], y[,2], use='pairwise.complete.obs',method=c('pearson'))
                }else{
                	cor=NA
                }
                c(rho= cor, rmse= rmse) })
            
            ## ============================== ##
            meanrhormse <- c()
            
            embed <- embed(c( rep(NA, 9+1), asv[libmat[r,1]:libmat[r,2]]), 10+1)
            meanpred <- embed[,-1]
              
            for(time.ahead in tp){
              
              n= sum(absasv[libmat[r,1]:libmat[r,2]]>0)
              y <- na.omit(cbind(asv[libmat[r,1]:libmat[r,2]], meanpred[, time.ahead]))
              rmse <-  sqrt(sum((y[,1]-y [,2])^2)/length(y)) / (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                
              meancor <- cor(y[,1], y[,2], use='pairwise.complete.obs',method=c('pearson'))
              meanrhormse <- rbind(meanrhormse, c(rho=meancor, rmse=rmse))
            }
            
            colname <- sprintf('%s_null_tp%s', colnames(ts)[k], tp)           
            dimnames(meanpred ) <- list(rowname, colname)
            ## ============================== ##
            
            stat <- rbind(stat, cbind(asv=colnames(ts)[k], replicate=r,tp=tp,
                                      AR=data.frame(t(rhormse)), mean=meanrhormse ))
            
            model <- rbind(model, cbind(ar=ar.pred, mean=meanpred))
        }
        
        list(stat, model)
    }    
   
    statics <- do.call(rbind, arres[,1])
    #plot(statics[statics[,2]==3,3:4])
    armodel <- c()
    for(j in tp){#j=1
    	
    	modeltmp <- sapply(arres[,2], function(x){ x[,j] })
    	modeltmp.names <- sapply(arres[,2], function(x){ colnames(x)[j] })
    	colnames(modeltmp) <- modeltmp.names
    	
    	armodel <- rbind(armodel, cbind(tp=j, modeltmp))
    	
    }
    
    meanmodel <- c()
    for(j in tp+10){#j=1
      
      modeltmp <- sapply(arres[,2], function(x){  x[,j] })
      modeltmp.names <- sapply(arres[,2], function(x){ colnames(x)[j] })
      colnames(modeltmp) <- modeltmp.names
      
      meanmodel <- rbind(meanmodel, cbind(tp=j-10, modeltmp))
      
    }
    ## ============================== ##
    
    predlist[[i]] <- armodel
    statall[[i]] <- statics
    meanmodels[[i]] <- meanmodel
    ## ============================== ##
}

saveRDS(meanmodels, 'Table/02_08_meanModel.rds')      
saveRDS(predlist, 'Table/02_08_ARModel.rds')    
saveRDS(statall, 'Table/02_08_AR_statics.rds')    

