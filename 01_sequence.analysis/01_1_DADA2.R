############################################################################
####
#### R script for Fujita (2019)
####
#### Bacteria sequence analysis by DADA2
#### Basically, 2018.4.18 Ushio created
#### 2020.10.31 Fujita edited
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('..')
#### setwd('~/Desktop/Microbiome_TimeSeries/MTS_5')
####
#### 2020.11.22. Changed refference
############################################################################
rm(list=ls())
# -- Record all processing
s <- Sys.time()

# -- Create directory to save
wd.dir <- '01_sequence.analysis'
save.dir <- sprintf('%s/01_output', wd.dir)
fastq.dir <- sprintf('%s/filtered_Fastq', wd.dir)
save.table <- sprintf('Table')

dir.create(save.dir)
dir.create(fastq.dir)
dir.create(save.table)
dir.create(sprintf('%s/01_4_AssignTaxnomy_claident', wd.dir))
dir.create(sprintf('%s/01_5_PICRUSt2', wd.dir))

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load library and function
library(dada2)
library(seqinr)
library(stringr)
library(ShortRead)
library(Biostrings)
library(ggplot2)

library.comment <- 
(sprintf('DADA2 v.%s
seqinr v.%s
stringr v.%s
ShortRead v.%s
Biostrings v.%s
ggplot2 v.%s\n',packageVersion("dada2"),packageVersion("seqinr"),packageVersion("stringr"),packageVersion("ShortRead"),packageVersion("Biostrings"),packageVersion("ggplot2")))

source( '01_sequence.analysis/function_refrence/MTS_function.R')

# -- Defined the name of forword and reverse sequence
pattern.f <- "__515f.forward.fastq.gz"
pattern.r <- "__515f.reverse.fastq.gz"

# -- Defined the name of reference data
reference1 <- "function_refrence/silva_nr99_v138_train_set_append_stdDNA.fa"
reference2 <- "function_refrence/silva_species_assignment_v138_append_stdDNA.fa"

# -- Difined the directory containing the fastq files
path <- 'Fastq' 
root.path <- list.files(path)
cutadpatpath="/home/your_user_name/.local/bin/cutadapt"
############################################################################

# -- Count remaining files
remain <- length(root.path)


# -- Processing DADA2 by each sequence run
for(f in root.path){ 
	object.start <- ls()
	# -- Record processing time
	cat(sprintf('Start inferring %s sequence \n',f)) ; time_start <- proc.time()[3]
	
	# -- Create directory to store filtered fastq files
	dir.create(sprintf('%s/%s', fastq.dir, f))
	pathR <- pathF <- paste(path,f, sep='/') 

	filtpathF <- file.path(sprintf('%s/%s', fastq.dir, f))
	filtpathR <- file.path(sprintf('%s/%s', fastq.dir, f))
	fastqFs <- sort(list.files(pathF, pattern = pattern.f ))
	fastqRs <- sort(list.files(pathR, pattern = pattern.r))

	if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

 	Check.primer(fwd= "NNNNNNGTGYCAGCMGCCGCGGTAA", rev ="NNNNNNGGACTACNVGGGTWTCTAAT", 
 		     path=pathF, patternf= pattern.f, patternr= pattern.r,
		     cutadpat=cutadpatpath )


	# -- Filtering low quality data
	cat('1. Filtering low quality samples...') ; filt_start <- proc.time()[3]
	out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
	              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
	              # if need remove Ns(3-6 bp) and  primers; 515F = 19 bp, 806rB = 20 bp
 	             # Output sequences of Claident already trimmed Ns and primers,
 	             truncLen=c(245,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
 	             compress=TRUE, verbose=TRUE, multithread=TRUE)
	cat(sprintf('Done ; elapsed %.02f\n', proc.time()[3]-filt_start))
	
	sample <- order(out[,1]-out[,2])[c(1:3,-1:1+round(median(1:nrow(out))),-2:0+nrow(out))]

	pdf(sprintf('%s/QualityProfile_%s_f.pdf', save.dir,f))
  	print(plotQualityProfile(file.path(pathF, fastqFs)[sample]))
  	dev.off()
	pdf(sprintf('%s/QualityProfile_%s_r.pdf', save.dir, f))
  	print(plotQualityProfile(file.path(pathR, fastqRs)[sample]))
  	dev.off()

	# -- Set up to infer Sequence Variants
	filtFs <- list.files(filtpathF, pattern= pattern.f, full.names = TRUE)
	filtRs <- list.files(filtpathR, pattern= pattern.r, full.names = TRUE)
	sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 4) # Assumes filename = samplename_XXX.fastq.gz
	sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 4) # Assumes filename = samplename_XXX.fastq.gz
	if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

	# --  Learn error rates
	cat('2. Learning error rates ...\n') ; er_start <- proc.time()[3]
	errF <- learnErrors(filtFs, nbases=1e10, multithread=TRUE)
	errR <- learnErrors(filtRs, nbases=1e10, multithread=TRUE)
	cat(sprintf('...Done ; elapsed %.02f\n', proc.time()[3]-er_start))
	
	pdf(sprintf('%s/plotErrors_%s_forward.pdf', save.dir,f))
  	print(plotErrors(errF, nominalQ=TRUE))
 	dev.off()
	pdf(sprintf('%s/plotErrors_%s_reverse.pdf', save.dir,f))
  	print(plotErrors(errR, nominalQ=TRUE))
 	dev.off()

	# --  Sample inference and merger of paired-end reads
	cat('3. inferring Sample and merge paired-end reads...') ; asv_start <- proc.time()[3]
	mergers <- vector("list", length(sample.names))
	names(mergers) <- sample.names
	for(sam in sample.names) {
    	derepF <- derepFastq(filtFs[[sam]])
    	ddF <- dada(derepF, err=errF, multithread=TRUE)
    	derepR <- derepFastq(filtRs[[sam]])
    	ddR <- dada(derepR, err=errR, multithread=TRUE)
    	merger <- mergePairs(ddF, derepF, ddR, derepR)
    	mergers[[sam]] <- merger
	}
	cat(sprintf('Done ; elapsed %.02f\n', proc.time()[3]-asv_start))

	# --  Construct sequence table and remove chimeras
	seqtab <- makeSequenceTable(mergers)
	saveRDS(seqtab, sprintf("%s/%s.rds", save.dir, f))
	
	remain <- remain-1
	cat(sprintf('Finish inferring %s and remain %s run\nelapsed time %.02f\n', f, remain, proc.time()[3]-time_start))
	object.end <- ls()
	rm(list=object.end[-which(object.end%in%object.start)])
	gc(reset=TRUE)
	gc(reset=TRUE)
}

############################################################################

gc(reset=TRUE)
gc(reset=TRUE)

# -- Merge multiple runs 
rds.file <- list.files(save.dir)
rds.file <- rds.file[grep(".rds",rds.file)]
for(i in 1:length(rds.file)) assign( letters[i], readRDS( sprintf('%s/%s', save.dir, rds.file[i]) ) )
inoculum <-  readRDS("00_original.inoculum/01_output/seqtab_pro.rds")
rownames(inoculum) <- paste(rownames(inoculum),"inoculum",sep="_")
st.all <- mergeSequenceTables(a,b,c,d,e,f,g)

# -- Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
saveRDS(seqtab, sprintf("%s/seqtab_all_with_inoculum.rds", save.dir))


# -- Assign taxonomy
cat('Start 1. Assign taxonomy and 2. addSpecies program...') ; asv_start <- proc.time()[3]
taxa <- assignTaxonomy(seqtab, sprintf("%s/%s", wd.dir,reference1), multithread=FALSE)
cat(sprintf('1. Done ; elapsed time %.02f\n', proc.time()[3]-asv_start))
taxa <- addSpecies(taxa[,-ncol(taxa)], sprintf("%s/%s", wd.dir,reference2))
cat(sprintf('2. Done ; elapsed time %.02f\n', proc.time()[3]-asv_start))

taxa.print <- taxa
head(taxa.print, n=10)
# --Excluding contamination sequence
seqtab.c <- seqtab

# -- Excluding Chloroplast and Mitochondria
taxa.print2 <- taxa.print[taxa.print[,'Kingdom']%in%c('Bacteria','Archaea'),]
mito = which(taxa.print2 == "Mitochondria", arr.ind=T)
chlo = which(taxa.print2 == "Chloroplast", arr.ind=T)
contami = which(taxa.print2 == "Escherichia-Shigella", arr.ind=T)

apply(taxa.print2, 2, table)
if(length(c(mito, chlo, contami)) > 0) taxa.print2 <- taxa.print2[-c(mito[,1], chlo[,1], contami[,1]),]
taxa.print2[is.na(taxa.print2)] <- 'unidentified'

# --Excluding standard DNA ASV
row.std <- grep('STD_', taxa.print2[,'Phylum'])
taxa.print.no.std <- taxa.print2[-row.std, ]
taxa.print.named <- as.data.frame(cbind(ID=sprintf("X_%04d", 1:nrow(taxa.print.no.std)), taxa.print.no.std))

# --Output ASV table
asv.ex.contami <- taxa.print.named
rownames(asv.ex.contami) <- asv.ex.contami[,1]
write.table(asv.ex.contami, file=sprintf("%s/Taxa_list.txt", save.table), sep="\t", quote=F, row.names=F)
saveRDS(asv.ex.contami,sprintf("%s/Taxa_list.rds", save.table))
saveRDS(cbind(asv.ex.contami,rownames(asv.ex.contami)),sprintf("%s/Taxa_list_with_sequence.rds", save.table))

# -- Output Fasta file
taxa.print.rename3 <- as.data.frame(asv.ex.contami)
fasta.table  <- cbind(as.character(taxa.print.rename3 $ID), rownames(taxa.print.named ))
	
write.fasta(as.list(fasta.table[,2]), fasta.table[,1],sprintf('%s/01_4_AssignTaxnomy_claident/01_Fasta_16S.fasta',wd.dir))
write.fasta(as.list(fasta.table[,2]), fasta.table[,1],sprintf('%s/01_5_PICRUSt2/01_Fasta_16S.fasta',wd.dir))

## ------------------------------------------------------------------ ##
# -- Out the sequence reads matrix
seqtab.bacteria <- seqtab[, rownames(taxa.print.named)]
if(all(colnames(seqtab.bacteria)==rownames(taxa.print.named))) {
	    colnames(seqtab.bacteria)= taxa.print.named[,1]
	}else{
	    stop(cat('colnames(seqtab.bacteria ) and rownames(asv.ex.contami) is different'))
}

rownames(seqtab.bacteria) <- paste('S_',rownames(seqtab.bacteria), sep='')
seqtab.bacteria2 <- cbind(ID=rownames(t(seqtab.bacteria)), t(seqtab.bacteria))
saveRDS(seqtab.bacteria, sprintf("%s/01_1_sequence_reads_Matrix.rds", save.table))
write.table(seqtab.bacteria2, sprintf('%s/01_5_PICRUSt2/Matrix.txt',wd.dir), sep="\t", quote=F, row.names=F)

# -- Out the matrix to check std reads proportion
seqtab.check <- seqtab
colnames(seqtab.check) <- taxa.print[,'Phylum']
	
sum.std <- rowSums(seqtab.check[,grep('STD_',colnames(seqtab.check))])
sum.na <-  rowSums(seqtab.check[,which(is.na(colnames(seqtab.check)))])
sum.pro <- rowSums(seqtab.check[,-c(grep('STD_',colnames(seqtab.check)),which(is.na(colnames(seqtab.check))))])
sum.reads <- rowSums(seqtab.check)
seqtab.check2 <- cbind(paste("S",rownames(seqtab.check),sep="_"),sum.reads,sum.std, sum.pro, sum.na, seqtab.check)
	
write.table(seqtab.check2, sprintf("%s/01_1_sequence_reads_check.txt", save.table), sep="\t", quote=F, row.names=F)

lf <- as.data.frame(rbind( cbind("Sum of STD DNA reads",sum.std), cbind("Sum of ASV reads",sum.pro),
		           cbind("Sum of NA sequence reads",sum.na),cbind("Sum of each sample reads",sum.reads)))
colnames(lf) = c("key","value")
gg <- ggplot(data=lf,aes(x=as.numeric(as.vector(value)))) +
		  geom_histogram(fill='white', color='indianred4',size=1)+
		  geom_hline(yintercept=0)+
		  facet_wrap(~key)+
		  theme_minimal(base_size=20)
		  
 	pdf(sprintf("%s/Reads.count.pdf",save.dir))
	plot(gg)	
	dev.off()

## ------------------------------------------------------------------ ##
# -- Save result

save.image(sprintf('%s/01_1_output.RData', save.dir, substr(Sys.time(), 1, 10)))

taxa.print3 <- taxa
rownames(taxa.print3) <- NULL
saveRDS(list(asv.table= taxa.print2, seqtab=seqtab[,rownames(taxa.print2)]), 'Table/for_STD_check.rds')
############################################################################
