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
source("functions/BidirectBlocklnlp_v0.7.4.R")

load.lib( c('rEDM', 'ggplot2', 'tidyr', 'cowplot'))

# -- Create directory to save
dir <- make.dir('06_community_prediction/06_02_output')

# -- Load data table
dlist  <- readRDS("Table/03_matrixList.rds")
libmat=dlist[[7]]
sml <- readRDS('Table/03_03_sample_info.rds')

# -- Decided regularized parameter range
E.range=1:20
theta_test <-  c(0.001, 0.01,0.05,0.1,0.2,0.5,1,2,4,8)
tp=1:10

for.parallel(8)

############################################################################

predlistT <- predlist0 <- statall <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ============================== ##
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ]

    ## ================================ ##
    smapres <- foreach(k = 1:ncol(ts), .combine=rbind)%dopar%{ #k=35
        
        asv <- ts[,k]
        absasv <- ts[,k]
            
        stat <- model <- c()
        for(r in 1:8){#r=1
  
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## -- Optimaization parameter
            E <- bestE_bidirect(asv, E_range = E.range, lib = libmat[-r,], criteria = "rmse", show_fig = F)
            
            ts_tmp <- embed(c(rep(NA,  E-1), asv),  E)
            theta <- bestT_bidirect_block_lnlp(ts_tmp, theta_range=theta_test, lib = libmat[-r,], criteria = "rmse", show_fig = F)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
            
            smapout <- s_map(asv, lib = libmat[-r,], pred=libmat[r,], 
                             E = E, theta= c(0, theta), tp = tp, tau = 1,
                             stats_only=FALSE)
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
            rhormse <- sapply(smapout$model_output, function(x){ #x=simplex.out$model_output[[1]]
                
                n= sum(absasv[libmat[r,1]:libmat[r,2]]>0)
                y <- na.omit(x[,2:3])
                rmse <-  sqrt(sum((y[,1]-y [,2])^2)/length(y)) / (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                
                if(!all(x[,2]==0, na.rm=TRUE)){
                  cor=cor(x[,2], x[,3], use='pairwise.complete.obs',method=c('pearson'))
                }else{
                  cor=NA
                }
                c(rho= cor, rmse= rmse) })
            
            stat <- rbind(stat, cbind(asv=colnames(ts)[k], replicate=r, smapout[,c(3, 1, 5)], t(rhormse)))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_smapT_tp%s', colnames(ts)[k], tp)
            
            tmpmatT <- matrix(NA, ncol=10, nrow=110, dimnames=list(libmat[r,1]:libmat[r,2], colname))
            rows=which(smapout$theta!=0)
            for(l in 1:length(rows)){ #l=1
                tmp <- smapout$model_output[[rows[l]]]
                rownames(tmp)[which(!is.na(tmp$time))] <- tmp$time[which(!is.na(tmp$time))]
                tmpmatT[,l] <- tmp[rownames(tmpmatT),'pred']
            }
            
            
            colname <- sprintf('%s_smap0_tp%s', colnames(ts)[k], tp)
            tmpmat0 <- matrix(NA, ncol=10, nrow=110, dimnames=list(libmat[r,1]:libmat[r,2], colname))
            rows=which(smapout$theta==0)
            for(l in 1:length(rows)){ #l=1
                tmp <- smapout$model_output[[rows[l]]]
                rownames(tmp)[which(!is.na(tmp$time))] <- tmp$time[which(!is.na(tmp$time))]
                tmpmat0[,l] <- tmp[rownames(tmpmat0),'pred']
            }
            
            rownames(tmpmatT) <- rownames(tmpmat0) <- rowname
            model <- rbind(model, cbind(tmpmatT, tmpmat0))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
        }
        
        list(stat, model[,grep('smapT', colnames(model))], model[,grep('smap0', colnames(model))] )
    }    
    
    ## ============================== ##
    
    statics <- do.call(rbind, smapres[,1])
    modelsmapT <- modelsmap0 <- c()
    for(j in tp){#j=1
      
        modelT <- sapply(smapres[,2], function(x){  x[,j]  })
        modelT.names <- sapply(smapres[,2], function(x){ colnames(x)[j] })
        colnames(modelT) <- modelT.names
        
        model0 <- sapply(smapres[,3], function(x){ x[,j]  })
        model0.names <- sapply(smapres[,2], function(x){ colnames(x)[j] })
        colnames(model0) <- model0.names
        
        modelsmapT <- rbind(modelsmapT, cbind(tp=j, modelT))
        modelsmap0 <- rbind(modelsmap0, cbind(tp=j, model0))
    }
    
    ## ============================== ##
    
    predlistT[[i]] <- modelsmapT
    predlist0[[i]] <- modelsmap0
    statall[[i]] <- statics
    
    ## ============================== ##
}

saveRDS(predlistT, 'Table/06_02_smapT_predicted_community')    
saveRDS(predlist0, 'Table/06_02_smap0_predicted_community')  
saveRDS(statall, 'Table/06_02_smap_statics')    

