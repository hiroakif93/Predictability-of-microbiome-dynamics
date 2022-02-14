############################################################################
####
#### R script for Fujita (2019)
####
#### Bacteria sequence analysis by DADA2
#### Basically, 2018.4.18 Ushio created
#### 2020.10.31 Fujita edited
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('..')
####
#### 2020.11.22. Changed refference
############################################################################

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('dada2', 'seqinr', 'stringr', 'ShortRead', 'Biostrings', 'ggplot2'))

# -- Create directory to save
dir = make.dir('00_ReadQC_to_Annotation/DADA2_Filtering2Annotatioin')

###########################################################################

# -- Defined the name of forword and reverse sequence
pattern.f <- "__515f.forward.fastq.gz"
pattern.r <- "__515f.reverse.fastq.gz"

# -- Defined the name of reference data
reference1 <- "Table/reference/silva_nr99_v138_train_set_append_stdDNA.fa"

# -- Difined the directory containing the fastq files
path <- 'MTS_fastq' 
root.path <- list.files(path)
root.path <- root.path[grep('Day', root.path)]

############################################################################

# -- Processing DADA2 by each sequence run
for(f in root.path){  #f= root.path[1]

    ## ============================================= ##
    
    # -- Create directory to store filtered fastq files
    dir.create(sprintf('%s/%s', dir$rdsdir, f))
    pathday <- paste(path, f, sep='/') 
    
    filtpathF <- file.path(sprintf('%s/%s', dir$rdsdir, f))
    filtpathR <- file.path(sprintf('%s/%s', dir$rdsdir, f))
    fastqFs <- sort(list.files(pathday, pattern = pattern.f, full.names = TRUE))
    pathF <-  paste(path, f, fastqFs, sep='/')
    fastqRs <- sort(list.files(pathday, pattern = pattern.r, full.names = TRUE))
    pathR <- paste(path, f, fastqRs, sep='/')
    
    if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
    
    ## ============================================= ##

	# -- Filtering low quality data
	cat('1. Filtering low quality samples...') ; filt_start <- proc.time()[3]

	## -- Filtered fastq directory
	filtFs <- file.path(filtpathF, "filtered", basename(fastqFs)) 
	filtRs <- file.path(filtpathR, "filtered", basename(fastqRs))
	
	out <- filterAndTrim(fwd=fastqFs, filt=filtFs,
	                     rev=fastqRs, filt.rev=filtRs,
 	             truncLen=c(245,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
 	             compress=TRUE, verbose=TRUE, multithread=TRUE)
	
	cat(sprintf('Done ; elapsed %.02f\n', proc.time()[3]-filt_start))

	sample <- order(out[,1]-out[,2])[c(1:3,-1:1+round(median(1:nrow(out))),-2:0+nrow(out))]

	pdf(sprintf('%s/QualityProfile_%s_f.pdf', dir$figdir,f))
  	print(plotQualityProfile( fastqFs[sample]))
  	dev.off()
	pdf(sprintf('%s/QualityProfile_%s_r.pdf', dir$figdir, f))
  	print(plotQualityProfile(fastqRs[sample]))
  	dev.off()

  	## ============================================= ##
  	# -- Set up to infer Sequence Variants
  	more0 <- out[which(out[,2]>0),]
  	
  	filtFs.ex0 <- filtFs[basename(filtFs)%in%rownames(more0)]
  	filtRs.ex0 <- filtRs[basename(filtFs)%in%rownames(more0)]
  	
  	sample.names <- sapply(strsplit(basename(filtFs.ex0), "_"), `[`, 4) # Assumes filename = samplename_XXX.fastq.gz
  	sample.namesR <- sapply(strsplit(basename(filtRs.ex0), "_"), `[`, 4) # Assumes filename = samplename_XXX.fastq.gz
  	if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  	names(filtFs.ex0) <- sample.names
  	names(filtRs.ex0) <- sample.names
  	
  	# --  Learn error rates
  	cat('2. Learning error rates ...\n') ; er_start <- proc.time()[3]
  	errF <- learnErrors(filtFs.ex0, nbases=1e10, multithread=TRUE)
  	errR <- learnErrors(filtRs.ex0, nbases=1e10, multithread=TRUE)
  	cat(sprintf('...Done ; elapsed %.02f\n', proc.time()[3]-er_start))
  	
  	pdf(sprintf('%s/plotErrors_%s.pdf', dir$figdir,f))
  	print(plotErrors(errF, nominalQ=TRUE))
  	print(plotErrors(errR, nominalQ=TRUE))
  	dev.off()
    
  	## ============================================= ##
  	
	# --  Sample inference and merger of paired-end reads
	cat('3. inferring Sample and merge paired-end reads...') ; asv_start <- proc.time()[3]
	mergers <- vector("list", length(sample.names))
	names(mergers) <- sample.names
	for(sam in sample.names) {
    	derepF <- derepFastq(filtFs.ex0[[sam]])
    	ddF <- dada(derepF, err=errF, multithread=TRUE)
    	derepR <- derepFastq(filtRs.ex0[[sam]])
    	ddR <- dada(derepR, err=errR, multithread=TRUE)
    	merger <- mergePairs(ddF, derepF, ddR, derepR)
    	mergers[[sam]] <- merger
	}
	cat(sprintf('Done ; elapsed %.02f\n', proc.time()[3]-asv_start))

	## ============================================= ##
	
	# --  Construct sequence table and remove chimeras
	names(mergers) <- sample.names
	seqtab <- makeSequenceTable(mergers)
	saveRDS(seqtab, sprintf("%s/%s.rds", dir$rdsdir, f))
	save.image(sprintf("%s/%s.RData", dir$rdatadir, f))
	
	gc(reset=TRUE)
	gc(reset=TRUE)
}

############################################################################

gc(reset=TRUE)
gc(reset=TRUE)

# -- Merge multiple runs 
rds.file <- list.files(dir$rdsdir, full.names = TRUE)
rds.file <- rds.file[grep(".rds",rds.file)]
rdslist <- lapply(rds.file, readRDS)

st.all <- mergeSequenceTables(tables=rdslist, repeats = "sum")

# -- Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
saveRDS(seqtab, "Table/00_seqtab_rmchimera.rds")


# -- Assign taxonomy
taxa <- assignTaxonomy(seqtab, sprintf("%s", reference1), multithread=FALSE)
saveRDS(taxa, "Table/00_taxaprint.rds")

############################################################################


