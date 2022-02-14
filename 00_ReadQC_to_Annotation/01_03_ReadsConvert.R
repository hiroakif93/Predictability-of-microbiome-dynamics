############################################################################
####
#### R script for Fujita (2019)
####
#### reads conversion
#### 2018.4.18 Ushio created
#### 2019.10.31 Fujita modified
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('..')
####
############################################################################

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('seqinr', 'ggplot2'))

# -- Load data table
taxa.print <- readRDS("Table/00_taxaprint.rds")
seqtab <- readRDS("Table/00_seqtab_rmchimera.rds")
sml <- read.table("Table/MTS_SampleSheet.txt",header=T)

############################################################################

# -- Nnumber of copies of standard DNA
std.copy.n <- c( 0.1, 0.05, 0.02, 0.01,0.005)*(6.02*10^13/1000000)/1000 *3

# -- Extract standard sequeces
detected.std.name <- unique(taxa.print[which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro"), "Phylum"])

n.std.seq <- which(substr(taxa.print[,"Phylum"], 1, 7) == "STD_pro")
std.table <- seqtab[,n.std.seq]
std.taxa <- taxa.print[n.std.seq, "Phylum"]

# --  STD reads - copy number relationship
# --  Rename colnames
colnames(std.table) <- std.taxa
# --  Merge the same colnames
new.std.table <- data.frame(std_rank1 = MergeSTD(detected.std.name[1], std.data = std.table),
                            std_rank2 = MergeSTD(detected.std.name[2], std.data = std.table),
                            std_rank3 = MergeSTD(detected.std.name[3], std.data = std.table),
                            std_rank4 = MergeSTD(detected.std.name[4], std.data = std.table),
                            std_rank5 = MergeSTD(detected.std.name[5], std.data = std.table))

# -- Interpolate 0 into missing value
new.std.table[is.na(new.std.table)] <- 0

# --  Linear regression
adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared
lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1]
r2.summary <- apply(new.std.table, 1, adj.r.fun)
coef.summary <- apply(new.std.table, 1, lm.coef.fun)
new.seqtab <- as.data.frame(seqtab[,-n.std.seq])

############################################################################

# -- Exclude NC w/o STD samples
nc.wo.std <- 0

# -- Convert sequence reads to copy numbers
coef.summary.all <- apply(new.std.table, 1, lm.coef.fun)
zero.coef <- which(!coef.summary.all > 0)
if(length(zero.coef) > 0) nc.wo.std <- c(nc.wo.std, zero.coef) # Add zero coefficient if any
coef.no.nc <- coef.summary.all[-nc.wo.std]
new.seqtab.no.nc <- as.matrix(new.seqtab[-nc.wo.std,])

# -- Coefficient check
all(coef.no.nc > 0)

# -- Duplicate check
rle (rownames(new.seqtab.no.nc)==names(coef.no.nc))

# -- Conversion
seqtab.conv0 <- new.seqtab.no.nc/coef.no.nc

# -- Add NC samples (with raw reads)
seqtab.invalid <- new.seqtab[nc.wo.std,]

# -- Combine valid & invalid samples
seqtab.conv.wo.std0 <- new.seqtab
seqtab.conv.wo.std0[c(1:nrow(new.seqtab))[-nc.wo.std],] <- seqtab.conv0
seqtab.conv.wo.std0[c(1:nrow(new.seqtab))[nc.wo.std],] <- seqtab.invalid

seqtab.conv <- seqtab.conv.wo.std0
taxa.wo.std <- taxa.print[-n.std.seq,]
dim(seqtab.conv)
dim(taxa.wo.std)

# -- Conversion of sample reads to calculated copy numbers
amp.factor <- 1/min(seqtab.conv[seqtab.conv != 0])
seqtab.conv.correct <- round(seqtab.conv*amp.factor) # Converted to integers

# -- Filt low quality sample (low quality means no-correlation with std.copy)
std.check <- apply(new.std.table, 1, function(x){  cor( x, as.matrix(std.copy.n))  })
low.q.sample <- which(std.check<0.7 | rowSums(new.seqtab) < 350)
length(low.q.sample)
seqtab.conv.correct.filt <- seqtab.conv.correct[-low.q.sample, ]
nrow(seqtab.conv.correct.filt)/10560
############################################################################

# --Excluding contamination sequence
seqtab.c <- seqtab.conv.correct.filt

# -- Excluding Chloroplast and Mitochondria
taxaBac <- taxa.print[taxa.print[,'Kingdom']%in%c('Bacteria','Archaea'),]
mito = which(taxaBac == "Mitochondria", arr.ind=T)
chlo = which(taxaBac == "Chloroplast", arr.ind=T)
contami = which(taxaBac == "Escherichia-Shigella", arr.ind=T) ## STD DNA vector

if(length(c(mito, chlo, contami)) > 0) taxaBac <- taxaBac[-c(mito[,1], chlo[,1], contami[,1]),]
taxaBac[is.na(taxaBac)] <- 'Unidentified'

# --Excluding standard DNA ASV
row.std <- grep('STD_', taxaBac[,'Phylum'])
taxa.print.no.std <- taxaBac[-row.std, ]
taxa.print.named <- as.data.frame(cbind(ID=sprintf("X_%04d", 1:nrow(taxa.print.no.std)), taxa.print.no.std))

# --Output ASV table
asv.ex.contami <- taxa.print.named
saveRDS(cbind(asv.ex.contami,rownames(asv.ex.contami)), "Table/Taxa_list_with_sequence.rds")

rownames(asv.ex.contami) <- NULL
rownames(asv.ex.contami) <- asv.ex.contami[,1]
saveRDS(asv.ex.contami, "Table/Taxa_list.rds")

# -- Output Fasta file
taxa.print.rename3 <- as.data.frame(asv.ex.contami)
fasta.table  <- cbind(as.character(taxa.print.rename3 $ID), rownames(taxa.print.named ))
	
write.fasta(as.list(fasta.table[,2]), fasta.table[,1],'Table/00_16S.fasta')

## ------------------------------------------------------------------ ##
# -- Out the sequence reads matrix
seqtab.bacteria <- seqtab.c[, rownames(taxa.print.named)]
if(all(colnames(seqtab.bacteria)==rownames(taxa.print.named))) {
	    colnames(seqtab.bacteria)= taxa.print.named[,1]
	}else{
	    stop(cat('colnames(seqtab.bacteria ) and rownames(asv.ex.contami) is different'))
}

rownames(seqtab.bacteria) <- paste('S_',rownames(seqtab.bacteria), sep='')
seqtab.bacteria2 <- cbind(ID=rownames(t(seqtab.bacteria)), t(seqtab.bacteria))
saveRDS(seqtab.bacteria, "Table/00_seqtabBac.rds")

# -- Out the matrix to check std reads proportion
seqtab.check <- seqtab
colnames(seqtab.check) <- taxa.print[,'Phylum']
	
sum.std <- rowSums(seqtab.check[,grep('STD_',colnames(seqtab.check))])
sum.na <-  rowSums(seqtab.check[,which(is.na(colnames(seqtab.check)))])
sum.pro <- rowSums(seqtab.check[,-c(grep('STD_',colnames(seqtab.check)),which(is.na(colnames(seqtab.check))))])
sum.reads <- rowSums(seqtab.check)
seqtab.check2 <- cbind(paste("S",rownames(seqtab.check),sep="_"),sum.reads,sum.std, sum.pro, sum.na, seqtab.check)
	
lf <- as.data.frame(rbind( cbind("Sum of STD DNA reads",sum.std), cbind("Sum of ASV reads",sum.pro),
		           cbind("Sum of NA sequence reads",sum.na),cbind("Sum of each sample reads",sum.reads)))
colnames(lf) = c("key","value")
gg <- ggplot(data=lf,aes(x=as.numeric(as.vector(value)))) +
		  geom_histogram(fill='white', color='indianred4',size=1)+
		  geom_hline(yintercept=0)+
		  facet_wrap(~key)+
		  theme_minimal(base_size=20)
		  
pdf(sprintf("%s/Reads.count.pdf",dir$figdir))
plot(gg)	
dev.off()

## ------------------------------------------------------------------ ##


############################################################################
