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
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[1]
    
    ## ============================== ##
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
    
    ## ============================== ##
    arres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=6

        asv <- tsscale[,k]
		absasv <- ts[,k]
        stat <- model <- c()
        
        for(r in 1:8){#r=1
            smlrep <- smlsub[smlsub$replicate.id==r,]
            ## ============================== ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_ar_tp%s', colnames(tsscale)[k], tp)
            
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
            tmp <- as.matrix(as.data.frame(ar.pred)[rownames(smlsub[smlsub$missing=='N', ]),])
            
            rhormse <- apply(tmp, 2, function(x){ #x=ar.pred[,1]
                
                nan <- which(smlrep$missing=='N')[which(which(smlrep$missing=='N')%in%c((max(20)+1):110))]
                x <- x[nan]
                
                n= sum(absasv[libmat[r,1]:libmat[r,2]][nan]>0)
                obs <- asv[libmat[r,1]:libmat[r,2]]
                y <- na.omit(cbind(obs[nan], x))
                
                rmse <-  sqrt(sum((y[,1]-y[,2])^2)/nrow(y)) # / mean(y[,1], na.rm=TRUE)
                rmse.weight <-  sqrt(sum((y[,1]-y [,2])^2)/n) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                rmse.weight <- ifelse(is.infinite(rmse.weight), NA, rmse.weight)
                
                if(!all(y[,1]==0, na.rm=TRUE)){
                	cor=cor(y[,1], y[,2], use='pairwise.complete.obs',method=c('spearman'))
                }else{
                	cor=NA
                }
                c(rho= cor, rmse= rmse, rmse.weight= rmse.weight) })
            
            ## ============================== ##
            meanrhormse <- c()
            
            statindex <- c(libmat[r,1]:libmat[r,2])[which(libmat[r,1]:libmat[r,2] %in% naindex)]
            embed <- embed(c( rep(NA, 10), asv[libmat[r,1]:libmat[r,2]]), 11)
            meanpred <- embed[,-1] #cbind(meanpred, rowMeans(embed[,-c(1:time.ahead)]))
            
            for(time.ahead in tp){ #time.ahead=1
              
              nan <- which(smlrep$missing=='N')[which(which(smlrep$missing=='N')%in%c((max(20)+1):110))]
              obs <- asv[libmat[r,1]:libmat[r,2]]
  
              n= sum(absasv[libmat[r,1]:libmat[r,2]][nan]>0)
              y <- cbind(obs, meanpred[, time.ahead])[nan,]
              rmse <-  sqrt(sum((y[,1]-y [,2])^2)/nrow(y))# / mean(y[,1], na.rm=TRUE)
              rmse.weight <-  sqrt(sum((y[,1]-y [,2])^2)/n) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
              
              meancor <- cor(y[,1], y[,2], use='pairwise.complete.obs',method=c('pearson'))
              meanrhormse <- rbind(meanrhormse, c(rho=meancor, rmse=rmse, rmse.weight= rmse.weight ))
            }
            
            colname <- sprintf('%s_null_tp%s', colnames(tsscale)[k], tp)           
            dimnames(meanpred) <- list(rowname, colname)
            ## ============================== ##
            
            stat <- rbind(stat, cbind(asv=colnames(tsscale)[k], replicate=r,tp=tp,
                                      AR=data.frame(t(rhormse)), null=meanrhormse ))
            
            model <- rbind(model, cbind(ar=ar.pred, mean=meanpred))
        }
        
        list(stat, model)
    }    
   
    statics <- do.call(rbind, arres[,1])
    #plot(statics[,4:5])
    armodel <- c()
    for(j in tp){#j=1
    	
    	modeltmp <- sapply(arres[,2], function(x){ #x = arres[,2][[1]]
    		x[,j]
    	})
    	armodel <- rbind(armodel, cbind(tp=j, modeltmp))
    	
    }
    
    meanmodel <- c()
    for(j in tp+10){#j=1
      
      modeltmp <- sapply(arres[,2], function(x){ #x = arres[,2][[1]]
        x[,j]
      })
      meanmodel <- rbind(meanmodel, cbind(tp=j-10, modeltmp))
      
    }
    ## ============================== ##
    
    predlist[[i]] <- armodel
    statall[[i]] <- statics
    meanmodels[[i]] <- meanmodel

    ## ============================== ##
}

saveRDS(meanmodels, 'Table/02_03_meanModel.rds')      
saveRDS(predlist, 'Table/02_03_ARModel.rds')    
saveRDS(statall, 'Table/02_03_AR_statics.rds')    
