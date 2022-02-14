############################################################################
####
#### R script for Fujita (2019)
####
#### Multivariate S-map predict all point
#### 2020.05 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('../')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
source("functions/BidirectSimplex_v0.7.4.R")
source("functions/BidirectBlocklnlp_v0.7.4.R")

load.lib( c('rEDM', 'ggplot2', 'tidyr', 'cowplot'))

# -- Create directory to save
dir <- make.dir('02_Predicting_microbiome_dynamics/Library_size')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
libmat=dlist[[7]]
sml <- readRDS('Table/sample_info.rds')
color <- readRDS('Table/color_palette.rds')

# -- Decided regularized parameter range
E.range <- 1:20
time.horizon <- 1:10
tp=c(1:10)
theta_test <-  c(0.001, 0.01,0.05,0.1,0.2,0.5,1,2,4,8)

for.parallel(detectCores())
############################################################################
statics.all <- c()

for(i in names(dlist)[-7]){ # for test : i = names(dlist)[6]
    
    cat( paste(i, '\n'))
    start <- Sys.time()
    ## ======================================================= ##
    ## -- Extracting one treatment matrix
    
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i,]
    
    ## ======================================================= ##
    ## -- Running in differnt library size
    statics.rep <- c()
    for(l in 1:7){# number of replicate as library : l=2
        
        cat( paste('library size', l, '...'))
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        ## -- Library replicates
        conb <- t(combn(1:8, l))
        lis <- c();  for(k in colnames(ts)){ lis <- rbind(lis, cbind(k, conb))}
        
        ## ======================================================= ##
        ## -- Estimating each ASV population dynamics
        
        tmp.result <- foreach(k = 1:nrow(lis), .combine=rbind, .packages=c("rEDM", "glmnet"))%dopar%{ # k=74
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            asv <- lis[k, 1]
            lib <- as.numeric(lis[k, -1])
            
            asv.scale <- as.numeric(scale(ts[,asv]))
            absasv <- ts[,asv]
            
            ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ## -- Running in each target replicate
            store.stat <- c()
            for(r in 1:8){ # r=3
                
                smlrep <- smlsub[smlsub$replicate.id==r,]
                ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## EDM parameters
                E <- bestE_bidirect(asv.scale, E_range = E.range, lib = libmat[lib,], criteria = "rmse", show_fig = F)		
                
                ts_tmp <- embed(c(rep(NA,  E-1), asv.scale),  E)
                theta <- bestT_bidirect_block_lnlp(ts_tmp, lib = libmat[lib,], criteria = "rmse", theta_range=theta_test, show_fig = F)
                
                ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## -- Simplex
                simplex.out <- simplex(asv.scale, lib = libmat[lib,], pred=libmat[r,], 
                                       E = E,  tp = time.horizon, tau = 1,
                                       stats_only=FALSE)
                
                rhormse <- sapply(simplex.out$model_output, function(x){ #x=simplex.out$model_output[[1]]
                    
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
                
                simpstat <- cbind(asv=asv, replicate=r, libsize=l, E=E, theta=theta,tp=1:10, simplex=as.data.frame(t(rhormse)))
                
                ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## -- S-map	
                smapout <- s_map(asv.scale, lib = libmat[lib,], pred=libmat[r,], 
                                 E = E, theta= c(0, theta), tp = tp, tau = 1,
                                 stats_only=FALSE)	       
                
                rhormseT <- sapply(smapout[smapout$theta!=0,]$model_output, function(x){ #x=simplex.out$model_output[[1]]
                    
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
                
                rhormse0 <- sapply(smapout[smapout$theta==0,]$model_output, function(x){ #x=simplex.out$model_output[[1]]
                    
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
                
                index=libmat[r, 1]:libmat[r,2]
                nan <- which(smlrep$missing=='N')[which(which(smlrep$missing=='N')%in%c((max(E.range)+1):110))]

                store.stat <- rbind(store.stat, cbind(treat=i, target.sample=sum(absasv[index[nan]]>0), simpstat,smapT=as.data.frame(t(rhormseT)), smap0=as.data.frame(t(rhormse0))))
                
            }#for end
            
            return(store.stat)
        } #foreach end
        
        statics.rep <- rbind(statics.rep, tmp.result)
    }
    
    saveRDS(statics.rep, sprintf("%s/libsize_dependency_%s.rds", dir$rdsdir, gsub("/", "_", i)))
    statics.all <- rbind(statics.all, statics.rep)
    
    ## ++++++++++++++++++++++++++++++++++++ ##
    end <- Sys.time()
    cat(sprintf('\nFinish.\n'))
    print(end-start)
}    

saveRDS(statics.all, "Table/libsize_dependency.rds")