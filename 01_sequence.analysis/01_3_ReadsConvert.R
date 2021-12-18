############################################################################
####
#### R script for Fujita (2019)
####
#### reads conversion
#### Basically, 2018.4.18 Ushio created
#### 2019.10.31 Fujita edited
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('..')
#### setwd('~/Desktop/Microbiome_TimeSeries/MTS_5')
####
############################################################################

# -- Create directory to save
dir.create("Table")
# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load workspace
load('01_sequence.analysis/01_output/01_1_output.RData')
load('01_sequence.analysis/01_output/01_2_output.RData')

# -- Delete large and unnecessary objects
rm(Fig.std) ;  rm(slope.summary) ; rm(new.seqtab2)
source('01_sequence.analysis/function_refrence/F1_HelperFunctions.R')
info <- read.table("Table/MTS_SampleSheet.txt",header=T)

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


############################################################################

# -- Save abundance matrix
mat <- seqtab.conv.correct.filt[,which(colnames(seqtab.conv.correct.filt) %in% rownames(taxa.print.named))]
if( all(colnames(mat) == rownames(taxa.print.named)) ) {colnames(mat) <- asv.ex.contami[,'ID']}else{stop('names is wrong')}
rownames(mat) <- paste("S",sep="_",rownames(mat))

mat2 <- cbind('sample'=rownames(mat), mat)
mat2 <- mat2[order(mat2 $sample),]

plot(mat[,10], type='l')

saveRDS(mat,"Table/Matrix_ab_01.rds")
write.table(mat2,sprintf('%s/Matrix_ab_01.txt',save.dir),sep="\t",quote= FALSE,row.names= FALSE)

# -- Save workspace

save.image(sprintf('%s/01_3_output.RData',save.dir))

############################################################################

table <- data.frame(row.names=c('Sequence reads (total)','Sequence reads (average)', 'ASVs'))
table$raw <- c(sum(st.all), mean(st.all[st.all>0]), ncol(st.all))
table[,'Removed chimera'] <- c(sum(seqtab), mean(seqtab[seqtab>0]), ncol(seqtab))
table[,'Bacteria sequence'] <- c(sum(seqtab.bacteria), mean(seqtab.bacteria[seqtab.bacteria>0]), ncol(seqtab.bacteria))
table[,'Standard DNA'] <- c(sum(new.std.table), mean(as.matrix(new.std.table)[as.matrix(new.std.table)>0]), ncol(new.std.table))

write.csv(table, 'Table/sequence result summary.csv', quote=F)

saveRDS(cbind(info, std_coef= coef.summary.all[ gsub('S_', '',as.character(info$sample))], std_cor= std.check[ gsub('S_', '',as.character(info$sample))]), 'Table/sample_info_01_3.rds')