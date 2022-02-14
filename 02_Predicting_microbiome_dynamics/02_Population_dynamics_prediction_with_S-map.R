############################################################################
####
#### R script for Fujita (2019)
####
#### Finding expected early warnings signal by visualizing
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
source("functions/BidirectSimplex_v0.7.4.R")
source("functions/BidirectBlocklnlp_v0.7.4.R")

load.lib( c('rEDM', "doParallel"))

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')

# -- Decided regularized parameter range
E.range=1:20
theta_test <-  c(0.001, 0.01,0.05,0.1,0.2,0.5,1,2,4,8)
tp=1:10

#for.parallel(detectCores())
cluster = makeCluster(detectCores()) 
registerDoParallel(cluster)
on.exit(stopCluster(cluster))
############################################################################

predlistT <- predlist0 <- statall <- c()
for(i in names(dlist)[-length(dlist)]){ # for test : i = names(dlist)[2]
    
    ## ++++++++++++ ##
    s <- proc.time()[3]
    ## ============================================= ##
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
    smapres <- foreach(k = 1:ncol(ts), .combine=rbind, .packages='rEDM')%dopar%{ #k=13
        
        asv <- tsscale[,k]
        absasv <- ts[,k]
            
        stat <- model <- c()
        for(r in 1:8){#r=5
            
            smlrep <- smlsub[smlsub$replicate.id==r,]
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
                
                nan <- which(smlrep$missing=='N')[which(which(smlrep$missing=='N')%in%c((max(E.range)+1):110))]
                index=libmat[r, 1]:libmat[r,2]
                x <- x[x$time%in% index[nan], ]
                n= sum(absasv[libmat[r,1]:libmat[r,2]][nan]>0)
                y <- x[,2:3]
                
                rmse <-  sqrt(sum((y[,1]-y [,2])^2)/nrow(y)) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                rmse.weight <-  sqrt(sum((y[,1]-y [,2])^2)/n) #/ (max(y[,1], na.rm=TRUE)-min(y[,1], na.rm=TRUE))
                rmse.weight <- ifelse(is.infinite(rmse.weight), NA, rmse.weight)
                
                if(!all(x[,2]==0, na.rm=TRUE)){
                    cor=cor(x[,2], x[,3], use='pairwise.complete.obs',method=c('spearman'))
                }else{
                    cor=NA
                }
                c(rho= cor, rmse= rmse, rmse.weight= rmse.weight) })
            
            stat <- rbind(stat, cbind(asv=colnames(tsscale)[k], replicate=r, smapout[,c(3, 1, 5)], t(rhormse)))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
            rowname <- rownames(smlsub[smlsub$replicate.id==r,])
            colname <- sprintf('%s_smapT_tp%s', colnames(tsscale)[k], tp)
            
            tmpmatT <- matrix(NA, ncol=10, nrow=110, dimnames=list(rowname , colname))
            rows=which(smapout$theta!=0)
            for(l in 1:length(rows)){ #l=1
                tmp <- smapout$model_output[[rows[l]]]
                tmpmatT[,l] <- tmp[,'pred']
            }
            plot(asv[445:554], type='l')
            lines(tmpmatT[,1], type='l')
            colname <- sprintf('%s_smap0_tp%s', colnames(tsscale)[k], tp)
            tmpmat0 <- matrix(NA, ncol=10, nrow=110, dimnames=list(rowname, colname))
            rows=which(smapout$theta==0)
            for(l in 1:length(rows)){ #l=1
                tmp <- smapout$model_output[[rows[l]]]
                tmpmat0[,l] <- tmp[,'pred']
            }
            
            model <- rbind(model, cbind(tmpmatT, tmpmat0))
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~ ##
        }
        
        list(stat, model[,grep('smapT', colnames(model))], model[,grep('smap0', colnames(model))] )
    }    
    
    ## ============================== ##
    
    statics <- do.call(rbind, smapres[,1])
    modelsmapT <- modelsmap0 <- c()
    for(j in tp){#j=1
        
        modelT <- sapply(smapres[,2], function(x){ #x = smapres[,2][[1]]
            x[,j]
        })
        modelT.names <- sapply(smapres[,2], function(x){ colnames(x)[j] })
        colnames(modelT) <- modelT.names
        
        model0 <- sapply(smapres[,3], function(x){ #x = smapres[,2][[1]]
          x[,j]
        })
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
    end <- proc.time()[3]
    cat( sprintf("%s %02f", i, end-s))
    ## ++++++++++++ ##
    
}

saveRDS(predlistT, 'Table/02_02_smapTModel.rds')    
saveRDS(predlist0, 'Table/02_02_smap0Model.rds')  
saveRDS(statall, 'Table/02_02_smap_statics.rds')    

## ============================== ##
