############################################################################
####
#### R script for Fujita (2019)
####
#### Shaping tables
#### 2020.9.28 Fujita
#### R 3.6.3
#### Set working directory of 'MTS' folder 
####  setwd('../')
############################################################################

# -- Create directory to save
dir.create('Table/for.mathematica')

# -- Load library, data table and function
source('functions/functions.R')
lib <- load.lib(c('vegan','tidyr'))

# -- Load data table
df <- readRDS('Table/00_seqtabBac.rds')
sml <- readRDS("Table/MTS_SampleSheet.rds")

############################################################################

## -- Threshold
replicate.th = 2
appearance.th = 10

############################################################################

dfnew <- smlnew <- c()
treat <- unique(sml[,2:3])
for(tr in 1:nrow(treat)){#tr=2
	
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
	smlsub <- Subset(data=sml, info=treat, row = tr) 
	smlsub$treat1 <- paste(treat[tr,2],treat[tr,1],sep='/')

	smlsub$treat2 <- paste(substr(treat[tr,2], 1, 1) ,
	                       substr(treat[tr,1], nchar(treat[tr,1]), nchar(treat[tr,1])),
	                       sep='')
	
	dfsub <- df[ rownames(df)%in%rownames(smlsub), ]
	
	## ========================================================== ##
	## -- Excluding low occurence species / pathway
	df.filt <- dfsub[,colSums(dfsub)>1]
    rowSums(df.filt)
	
    ## ========================================================== ##
	# -- Insert NA 
	smlsub <- smlsub[order(smlsub$replicate.id),]
	df.filt <- df.filt[rownames(smlsub),] 
    rownames(df.filt) <- rownames(smlsub)
	
    smlsub$missing <- 'N'
	smlsub[ which( is.na(df.filt[,1]) ), 'missing'] <- 'Y'

	# -- Interpotion
	for(r in 1:8){ #r=1
	    
	    smlrep <- smlsub[smlsub$replicate.id==r, ]
	    if(smlsub[1, 'missing']=='Y') df.filt[ rownames(smlsub)[1], ] <- df.filt[ rownames(smlsub)[2], ]
	    if(smlsub[110, 'missing']=='Y') df.filt[ rownames(smlsub)[110], ] <- df.filt[rownames(smlsub)[109], ]
	    
	}
	
	if( length( which( is.na(df.filt[,1]))) > 0 ){
	  df.int <- Interpotion(df.filt,method='mean') # method is 'mean' or 'linear'
	}else{
	  df.int <- df.filt
	}
	
	## ========================================================== ##
	## -- For Energy landscape matrix
	df.math <- as.data.frame(t(apply(df.filt,1,function(x){ x/sum(x)})))

	math 	<- cbind('time'=as.numeric( gsub('Day_', '', smlsub$time)), 'replicate'= smlsub$replicate.id, df.math)
	write.csv(math[order(math$time),],sprintf('Table/for.mathematica/ra_base_all.rep_%s.csv', unique( smlsub$treat2 )),row.names=FALSE)

	## ========================================================== ##
	## -- For EDM matrix, insert NAs among replicates
	
	## -- Filtering ASVs and statics
	ts.num <- smlsum <- c()		
	for(r in 1:8){ #r=4
	    
	    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
	    ## -- Extracting ASVs which appeared more than 'appearance.th'
	    smlrep <- smlsub[smlsub $replicate.id==r, ]
	    ts.rep <- df.int[rownames(smlrep), ]
	    bi <- ts.rep; bi[bi>0] <- 1
	    ts.num <- rbind(ts.num, colSums(bi)> appearance.th)
	}
	
	ts.filt <- df.int[ ,colSums(ts.num) >= replicate.th]
	    
	d.list <- lapply(1:8, function(x){ rbind(ts.filt[ rownames( smlsub[ which(smlsub$replicate.id==x), ] ),], NA) } )
	edmmat <- do.call(rbind, d.list)
	## ========================================================== ##	
	dfnew[[unique( smlsub$treat1 )]] <- edmmat
	smlnew <- rbind(smlnew,  smlsub) 
	
}

smlnew$timeChar <- smlnew$time
smlnew$time <- as.numeric( gsub('Day_', '', smlnew$time))
smlnew$treat1 <- factor(smlnew$treat1, levels=unique(smlnew$treat1))

saveRDS(dfnew,'Table/matrixList.rds')
saveRDS(smlnew,'Table/sample_info.rds')

############################################################################
