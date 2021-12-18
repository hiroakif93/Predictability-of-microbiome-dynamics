############################################################################
####
#### R script for Fujita (2019)
####
#### Shaping tables
#### 2020.9.28 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6')
####  setwd('..')
############################################################################

# -- Create directory to save
dir <- make.dir('02_shaping_data.table/02_output')
dir.create('Table/for.mathematica')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load library, data table and function
lib <- load.lib(c('vegan','tidyr','ggplot2'))
source('02_shaping_data.table/02_function.R')

# -- Load data table
df <- readRDS('Table/Matrix_ab_01.rds')
info <- readRDS("Table/sample_info_01_5.rds")
func <- as.data.frame(t(read.delim('01_sequence.analysis/01_5_PICRUSt2/output/pathways_out/path_abun_unstrat.tsv.gz',header=T,row.names=1)))

# -- The object to store all result
fn.al <-info.al<- df.al <- c()

############################################################################
# -- Run by treatment
treat <- unique(info[,3:5])
for(tr in 1:nrow(treat)){#tr=1
	
	# -- Subset one treatment data
	info2 		 <- Subset(data=info,info=treat, row = tr) 
	info2$treat1 <- paste(treat[tr,1],treat[tr,2],treat[tr,3],sep='_')
	
	t.name2 <-  c(); t.name <- substr( sapply(strsplit( unique( info2$treat1 ), '-','_'),'['), 1,1)
	for(i in 1:length(t.name)){ t.name2 <- paste(t.name2, t.name[i], sep='') }
	info2$treat2 <- paste(t.name2,substr(as.character(treat[tr,2]),1,1),substr(as.character(treat[tr,3]),1,1),sep='')
	
	df2 <- df[ rownames(df)%in%info2[,1], ]
	func2 <- func[rownames(func)%in%info2[,1], ]
	
	# -- Excluding low occurence species / pathway
	df.filt <- df2[,colSums(df2)>1]
	fn.filt <- func2[,colSums(func2)>1]

	# -- Insert NA 
	info3 <- info2[order(info2 $replicate.id),]
	df3   <- df.filt[as.character(info3[,1]),] 
	func3 <- fn.filt[as.character(info3[,1]),]
	
	rownames(func3) <- rownames(df3) <- info3[,1]

	info3 $missing <- 'N'
	info3[ which( is.na(df3[,1]) ), 'missing'] <- 'Y'
	
	# -- Counting species richness
	binary <- df3
	binary[binary > 0] <- 1
	
	info3$richness <- rowSums(binary)
	
	# -- Interpotion
	for(r in 1:8){ #r=1
	    
	    info3sub <- info3[info3$replicate.id==r, ]
	    if(info3sub[1, 'missing']=='Y') df3[ as.character(info3sub[1, 'sample']), ] <- df3[ as.character(info3sub[2, 'sample']), ]
	    if(info3sub[110, 'missing']=='Y') df3[ as.character(info3sub[110, 'sample']), ] <- df3[as.character(info3sub[109, 'sample']), ]
	    
	}
	
	if( length( which( is.na(df3[,1]))) > 0 ){
	  df.int <- Interpotion(df3,method='mean') # method is 'mean' or 'linear'
	}else{
	  df.int <- df3
	}
	
	if( length( which( is.na(func3[,1]))) > 0 ){
	  func.int <- Interpotion(func3,method='mean')
	}else{
	  func.int <- func3
	}
	df.math <- as.data.frame(t(apply(df3,1,function(x){ x/sum(x)})))
	fn.math <- as.data.frame(t(apply(func3,1,function(x){ x/sum(x)})))
	
	math 	<- cbind('time'=as.numeric(info3$time), 'replicate'= info3$replicate.id, df.math)
	write.csv(math[order(math$time),],sprintf('Table/for.mathematica/ra_base_all.rep_%s.csv', unique( info2$treat2 )),row.names=FALSE)
	
	math2 	<- cbind('time'=as.numeric(info3$time), 'replicate'= info3$replicate.id, fn.math)
	write.csv(math2[order(math2$time),],sprintf('Table/for.mathematica/function_base_all.rep_%s.csv', unique( info2$treat2 )),row.names=FALSE)
	
	# -- Insert NA 
	df.int.temp 		<- as.data.frame(t(df.int))
	df.int2 			<- t( df.int.temp[ colnames(df2) , ] )
	colnames(df.int2) 	<- colnames(df2)
	
	func.int.temp 		<- as.data.frame(t(func.int))
	func.int2 			<- t( func.int.temp[ colnames(func2) , ] )
	colnames(func.int2 ) 	<- colnames(func2)
	
	df.int2[is.na(df.int2)] <- 0
	func.int2[is.na(func.int2)] <- 0
	
	df.al <- rbind(df.al, df.int2) ; info.al <- rbind(info.al, info3) ; fn.al <- rbind(fn.al, func.int2)
	
}

############################################################################

# --Calculate alpha diversity index
shannon <- as.matrix(diversity(df.al,'shannon'))
info.al$asv.shannon <- shannon
simpson <- as.matrix(diversity(df.al,'simpson'))
info.al$asv.simpson <- simpson
shannon2 <- as.matrix(diversity(fn.al,'shannon'))
info.al$function.shannon <- shannon2
simpson2 <- as.matrix(diversity(fn.al,'simpson'))
info.al$function.simpson<- simpson2

# -- Check structure
no.na.df <- df.al[which(info.al $missing=='N'), ]
if( all(dim(no.na.df) == dim(df)) == FALSE  ) stop('dimention is wrong')
if( all(rownames(no.na.df)%in%rownames(df)) == FALSE ) stop('rownames is wrong')
df.al2 <- df.al[order(rownames(df.al)),]
fn.al2 <- df.al[order(rownames(fn.al)),]

# -- Change name
info.al$resource2 <- gsub('^Oat-Peptone$','Medium B',
						gsub('^Oat$','Medium A',
						gsub('^Peptone$','Medium C',
									info.al$resource )))
info.al$resource2 <- factor(info.al$resource2, levels=c('Medium A','Medium B','Medium C'))
info.al$treat1 <- gsub('^OSN','Soil/medium A',
					gsub('^OWN','Water/medium A',
					gsub('^OPSN','Soil/medium B',
					gsub('^OPWN','Water/medium B',
					gsub('^PSN','Soil/medium C',
					gsub('^PWN','Water/medium C', info.al$treat2 ))))))
info.al <- info.al[info.al$perturbation=='N',]

# -- Save data table and Rdata
saveRDS(df.al[order(rownames(df.al)),], 'Table/Matrix_int_02_1.rds')
saveRDS(fn.al[order(rownames(fn.al)),], 'Table/Function_int_02_1.rds')
saveRDS(info.al[order(info.al$sample),], 'Table/sample_info_02_1.rds')

# -- Save working space
save.image(sprintf('%s/02_data_shaping.RData',dir$rdatadir))

############################################################################
## ---------------------------------------------------------------------- ##
# -- Save environment
sink(file=sprintf("%s/session_01.txt",dir$dir))
print(Sys.time())
sessionInfo()
sink()
