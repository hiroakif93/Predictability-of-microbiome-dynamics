############################################################################
####
#### R script for Fujita (2019)
####
#### Color palette of each taxanomy level
#### 2020.9.28 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_5')
#### setwd('..')
############################################################################
dir <- make.dir('02_shaping_data.table/02_output')

# -- Load library and function
lib <- load.lib(c('RColorBrewer'))
source('02_shaping_data.table/02_function.R')

# -- Load data table and setting parameters
d		<- readRDS("Table/Matrix_int_02_1.rds")			# otu-abundance table
taxa 	<- as.data.frame(readRDS("Table/Taxa_list.rds"))	# color palette table
taxaLabel 	 		<- c("Phylum","Class","Order","Family","Genus", "identified")
taxa.palette 		<- vector('list',length=length(taxaLabel))
names(taxa.palette) <- taxaLabel

# -- Select color manually
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' ,]
qual_col_pals <- qual_col_pals[c('Set1','Dark2','Accent','Set2','Set3', 'Pastel1', 'Pastel2'),]
col.vector = c( brewer.pal(8, 'Set1'), brewer.pal(7, 'Dark2'), brewer.pal(7, 'Accent'),brewer.pal(7, 'Set2'),
				brewer.pal(8, 'Set3'), brewer.pal(8, 'Pastel1'), brewer.pal(7, 'Pastel2'))

unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

palette <- c(unique(col.vector),
			 'burlywood1','lavender','aliceblue')
#plot(1:length(palette), col= palette, pch=19, cex=4)
############################################################################

taxa$identified <- NA

taxa$identified <- 	apply(taxa, 1, function(x){ # x=taxa[4,]
		
						unident <- which(x=='unidentified')
						
						if(length(unident)>0){
							last.ident <- x[min(unident)-1]
						}else{
							last.ident <- x['Species']
						}
						
						last.ident	})

saveRDS(taxa, "Table/Taxa_list_02.rds")						
############################################################################

# -- Colored by taxanomy level
rownames(taxa) <- taxa[,1]

for(i in 1:length(taxaLabel)){
    
	df <- Taxa.mat(x=d, y=taxa, taxaLabel[i])

	# -- Defined low appearance ASV as 'Others'
	sort <- df[order(rowSums(df), decreasing=TRUE),]
	if(nrow(sort) > length(palette)){
		temp <- colSums(sort[ (length(palette)+1):nrow(sort) ,])
		sort <- rbind( sort[1:length(palette),], temp)
		rownames(sort) <- c(rownames(sort[1:length(palette),]), 'Others')
	}
	
	# -- Defined color, but Unidentified and Others is grey30 or grey60
	sort 		<- sort[rownames(sort) != 'unidentified' & rownames(sort) != 'Others',]
	taxa$color 	<- NA
	for( j in 1:length(rownames(sort)) )taxa[ taxa[,taxaLabel[i]] == rownames(sort)[j],]$color <- rev(palette[j])

	taxa$color[ which(taxa[, taxaLabel[i]] == 'unidentified') ] <- 'grey60'
	taxa$color[ is.na(taxa$color) ] <- 'grey30'
	
	color <- unique(taxa[, c( which(colnames(taxa)== taxaLabel[i]),ncol(taxa))]) 
	
	rownames(color) 	<- color[,1]
	taxa.palette[[i]] 	<- color
	
}


############################################################################

# Save palette

saveRDS(taxa.palette,'Table/taxa_color_02.rds')
saveRDS(taxa.palette,'02_shaping_data.table/02_output/taxa_color_02.rds')

############################################################################
sink(file=sprintf("%s/session_02.txt",dir$dir))
print(Sys.time())
sessionInfo()
sink()
