############################################################################
####
#### R script for Fujita (2019)
####
#### Compiling data matrix
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('RColorBrewer', 'circlize', 'ggplot2', 'tidyr', 'ggh4x', 'vegan'))
dir <- make.dir('03_dynamics_indices/01_compile')

## -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

## -- Load data table and setting parameters
ts <- as.data.frame(readRDS("Table/Matrix_int_02_1.rds"))		# ASV-abundance table
info <- readRDS("Table/sample_info_02_1.rds")[,1:13]
rownames(info) <- info[,1]
sml <- info[order(info$replicate.id),-1]
taxa <- readRDS('Table/Taxa_list_02.rds')

## -- Threshold
replicate.th = 1
appearance.th = 5

dir.create('Table/mathmatica')

############################################################################

dim(ts[rownames(info), colSums(ts[rownames(info), ]>0)>0])
write.csv(ts[rownames(info), colSums(ts[rownames(info), ]>0)>0], 'Table/No-perturbation_community_matrix.csv')

############################################################################

df.list <- asv.list <- list()
smllist <- asvname <- c()

pdf(sprintf('%s/ASV_frequency.pdf', dir$figdir), w=10, h=15)
for(i in as.character(unique(info$treat1))){ # for test : i=as.character(unique(info$treat1))[1]
	
	## ============================================================ ##
	## -- Extracting one treatment matrix
	smlsub <- sml[sml $treat1==i,]
	ts.sub <- ts[rownames(smlsub),]
	ts.sub <- ts.sub[,colSums(ts.sub)>0] 
	taxa.sub <- taxa[colnames(ts.sub), ]
	
	## ============================================================ ##
	## -- Filtering ASVs and statics
	ts.num <- smlsum <- c()		
	for(r in 1:8){ #r=4
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Extracting ASVs which appeared more than 'appearance.th'
		smlrep <- smlsub[smlsub $replicate.id==r, ]
		ts.rep <- ts.sub[rownames(smlrep), ]
		bi <- ts.rep; bi[bi>0] <- 1
		ts.num <- rbind(ts.num, colSums(bi)> appearance.th)
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Counting dominant ASV
		smlrep$richness <-  rowSums(bi)
		taxa.sub[colnames(bi),sprintf('num.sample_r%s', r)] <- colSums(bi)
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Mean cv
		asv.cvs <-  apply(ts.rep[, which(colSums(bi)> appearance.th)], 2, function(x){
						dly <- embed( c(rep(NA, 9), x), 10)
						mean <- rowMeans(dly); sd <- apply(dly, 1, sd)
						return( sd/mean )		}) 
		smlrep$abs.cvs <- rowMeans(asv.cvs, na.rm=TRUE)
		taxa.sub[colnames(ts.rep),sprintf('cv_r%s', r)] <- apply(ts.rep, 2, function(x) sd(x)/mean(x))
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Total abundance
		smlrep$total.abundacne <- rowSums(ts.rep[, which(colSums(bi)> appearance.th)])
		dly <- embed( c(rep(NA, 9), smlrep$total.abundacne), 10)
		smlrep$ta.cv <- apply(dly, 1, sd)/rowMeans(dly)
		taxa.sub[colnames(ts.rep),sprintf('frequency_r%s', r)] <- apply(ts.rep/rowSums(ts.rep), 2, max)
		
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		## -- Relative abundance cv
		rel <- ts.rep[, which(colSums(bi)> appearance.th)]/rowSums(ts.rep[, which(colSums(bi)> appearance.th)])
		## -- Mean cv
		asv.cvs <-  apply(rel, 2, function(x){
						dly <- embed( c(rep(NA, 9), x), 10)
						mean <- rowMeans(dly); sd <- apply(dly, 1, sd)
						return( sd/mean )		})
		smlrep$rel.cvs <- rowMeans(asv.cvs, na.rm=TRUE)
					
		## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
		
		smlsum <- rbind(smlsum , smlrep)
	}
	
	ts.filt <- ts.sub[ ,colSums(ts.num) >= replicate.th]
	print(dim(ts.sub)-dim(ts.filt))
	taxa.filt <- taxa.sub[colnames(ts.filt), ]
	taxa.filt[, 'num.replicate'] <- colSums(ts.num)[rownames(taxa.filt)]
	taxa.filt[colnames(ts.filt),'frequency_all'] <- apply(ts.filt/rowSums(ts.filt), 2, max, na.rm=TRUE)
	taxa.filt[colnames(ts.filt),'cv_mean_all'] <- rowMeans(taxa.filt[,grep('cv_', colnames(taxa.filt))], na.rm=TRUE)
	
	## ============================================================ ##
	## -- Taxonomy property
	tmp <- taxa.filt[,grep('frequency_r', colnames(taxa.filt))]
	asvorder <- apply(taxa.filt[,c('ID', 'Family')], 1, paste, collapse='-')[hclust(dist(tmp) )$order]
	rorder <- c(colnames(tmp)[hclust(dist(t(tmp)) )$order],'frequency_all')
	
	lf <- gather(cbind(id=apply(taxa.filt[,c('ID', 'Family')], 1, paste, collapse='-'), taxa.filt[,grep('frequency', colnames(taxa.filt))]), key, value, -1)
	lf$id <- factor(lf$id , level= asvorder)
	lf$key <- factor(lf$key , level=rorder)
	
	lf$group <- 'each'; lf$group[grep('_all', lf$key)] <- 'all'
	g <- ggplot(lf)+
		geom_tile(aes(x=key, y=id, fill=value), color='grey90')+
		scale_fill_gradientn(colors=brewer.pal(11,'Spectral'))+
		facet_wrap(~group, scales='free_x')+
		force_panelsizes(cols = c(1/8, 1))+
		labs(title=i)+
		theme_minimal()
	plot(g)
	## ============================================================ ##
	## -- For the EDM, NA inserting
	d.list <- lapply(1:8, function(x){ rbind(ts.filt[ rownames( smlsum[ which(smlsum$replicate.id==x), ] ),], NA) } )
	new.matrix <- do.call(rbind, d.list)
	asvname <- c(asvname, colnames(new.matrix))
	print(i); print(dim(new.matrix))
	write.csv(new.matrix, sprintf('Table/%s_community_matrix.csv', gsub('/', '_', i)))
	
	df.list[[i]] <- new.matrix
	smllist <- rbind(smllist, smlsum)
	asv.list[[i]] <- taxa.filt
	
	df <- as.data.frame(ts.filt)
	smltmp <- smlsum[order(smlsum$time), ]
	df[rownames(smltmp[smltmp $missing=='Y',]),] <- NA
	write.csv(cbind(time=as.numeric(smltmp[,'time']), replicate= smltmp[,'replicate.id'], df[rownames(smltmp),]/rowSums(df[rownames(smltmp),])),
				sprintf('Table/mathmatica/ocdata_%s.csv', gsub('/', '_', i)), row.names=FALSE)
}
dev.off()

## -- Store the result
na.row <- which(is.na(df.list[[1]][,1]))
lib_mat <- cbind(na.row-110,na.row-1)
df.list[['lib_mat']] <- lib_mat

saveRDS(df.list,'Table/03_matrixList.rds')
saveRDS(smllist,'Table/03_sample_info.rds')
saveRDS(asv.list,'Table/03_asv_info.rds')

############################################################################
if(FALSE){
mat <- matrix(0, ncol=length(unique(asvname)), nrow=5280, dimnames=list(rownames(sml), unique(asvname)))

for(i in names(df.list)[-7]){ # for test : i=as.character(unique(info$treat1))[1]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    tssub <- as.matrix(na.omit(df.list[[i]]))
    mat[rownames(tssub), colnames(tssub)] <- tssub
    
}    

nmds <- metaMDS(mat/rowSums(mat), distance='bray',
                parallel=5)

df <- cbind(sml, nmds$points[rownames(sml), ])
saveRDS(df, 'Table/03_01_NMDS_all_Sample.rds')
}

############################################################################

## -- Making color palette
taxa <- as.matrix(readRDS('Table/Taxa_list_02.rds'))
bin_info <- unique(read.table('Table/08_metagenome/bin_summary.txt', sep='\t', header=T)) ; dim(bin_info)

taxa['X_0003', 2:8] <- as.character(as.matrix(bin_info)[22,4:10])
## -- ASV Color

sep <- 5
col1 <- colorRamp2( c(1, ceiling(sep), sep*2), c('white', 'red', 'firebrick3') ) # plot( c(1: (sep*2))[-1]~1, col=col1(1: (sep*2))[-1] , pch=19)
col2 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'orange', 'darkorange2') )
col3 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'gold', 'darkgoldenrod1') )
col4 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'olivedrab2', 'chartreuse4') )
col5 <- colorRamp2( c(1, ceiling(sep/2), 6), c('white', 'steelblue1', 'skyblue3') )
col6 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'royalblue1', 'dodgerblue4') )
col7 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'plum3', 'slateblue2') )
col8 <- colorRamp2( c(1, ceiling(sep/2), sep), c('white', 'peachpuff', 'peachpuff4') )

pal <- c(col1(seq(1, (sep*2+1), 2))[-1],  col6(seq(1, (sep*2+1), 2))[-1], col5(seq(1, (sep*2+1), 2))[-1],  col4(seq(1, (sep*2+1), 2))[-1], col7(seq(1, (sep*2+1), 2))[-1], 
         col2(seq(1, (sep*2+1), 2))[-1], col3(seq(1, (sep*2+1), 2))[-1], col8(seq(1, (sep*2+1), 2))[-1] )
length(pal)
colmat <- matrix(pal, ncol=8, nrow= sep, byrow=TRUE); colmat <- t(apply(colmat, 1, rev))
colpal <- c( as.vector(colmat), brewer.pal(12, 'Set3')[-9] )
lf <- gather( as.data.frame( cbind(1:sep, sapply(1:7, function(x) {rep(x, sep)}) ) ), key,value, -1 )
plot(lf[,c(1,3)], col= colpal, pch=19) ;dev.off()

colmat <- rbind(brewer.pal(9, 'YlOrRd')[3:7], brewer.pal(9, 'YlGnBu')[3:7], brewer.pal(9, 'YlGn')[3:7], 
				brewer.pal(9, 'Purples')[3:7], brewer.pal(9, 'PuRd')[3:7], c('cornsilk1','cornsilk2', 'bisque2', 'bisque3', 'bisque4'))
colmat <- t(apply(colmat, 1, rev))
colpal <- c( as.vector(colmat)) 

lf <- gather( as.data.frame( cbind(1:5, sapply(1:5, function(x) {rep(x, sep)}) ) ), key,value, -1 )
plot(lf[,c(1,3)], col= colpal, pch=19) ;dev.off()
plot(1:length(colpal), col=colpal, pch=19)

colorlist <- c()
for(i in colnames(taxa)[1:7][-2]){ # i=colnames(taxa)[1:7][-2][7]
	
	taxaname <- as.character(unique(taxa[,i]))
	abundancemat <- as.data.frame(matrix(0, ncol=2, nrow=length(taxaname), dimnames=list(taxaname, NULL)))
	
	for(l in names(df.list)[-7]){ #l=names(df.list)[1]
		
		ts <- t(Taxa.mat(df.list[[l]], taxa, i))
		ts <- ts/rowSums(ts)
		dominant <- sort(colSums(ts, na.rm=TRUE))
		
		abundancemat[names(dominant), 1] <- abundancemat[names(dominant), 1] + dominant

	}

	abundancemat <- abundancemat[order(abundancemat[,1], decreasing=TRUE), ]
	abundancemat[,2] <- 'grey30'
	
	if( sum(rownames(na.omit(abundancemat)) == 'unidentified' ) > 0 ){
	    abundancemat[ 'unidentified' ,2] <- 'grey90'
	}
	if(i=='Genus') { abundancemat[taxa['X_0003', 7],2] <- 'yellow'}
	if(any(rownames(abundancemat)%in%taxa['X_0003', c(1,5:7)]))abundancemat[rownames(abundancemat)%in%taxa['X_0003', c(1,5:7)],2]<- 'yellow'

	abundancemat[-c(which(rownames(abundancemat)=='unidentified'), which(rownames(abundancemat)%in%taxa['X_0003', c(1,5:7)])),][1:length(colpal),2] <- colpal

	pal <- na.omit(abundancemat)[,2] ; names(pal) <- rownames(na.omit(abundancemat))
	colorlist[[i]] <- pal
	
	###########################################
	
	names(pal)[pal=='grey30'] <- 'Others'
	lf <- as.matrix(cbind(names=rownames(abundancemat), as.data.frame(abundancemat)))
	lf[lf[,3]=='grey30','names'] <- 'Others'
    lf <- as.data.frame(unique(lf))
    
    if(sum(unique(names(pal))=='unidentified')>0){	
        
        gcol <- pal[c(unique(names(pal))[-which(unique(names(pal))=='unidentified')],'unidentified')]
        lf$names <- factor(lf$names, levels=c(unique(names(pal))[-which(unique(names(pal))=='unidentified')],'unidentified') )
    }else{
        gcol <- pal[c(unique(names(pal)))]
    }

	g <- ggplot(lf)+
	    geom_point(aes(x=names, y=names, color=names), size=8)+
	    scale_color_manual(values=gcol)+
	    theme_bw(base_size=27)+
	    theme(legend.position='bottom')+
	    guides(color=guide_legend(nrow=5) )
	pdf(sprintf('%s/color_legend_%s.pdf', dir$figdir, i), w=35, h=15); plot(g); dev.off()
	###########################################
	
}

## -- replicate Color
rep.col <- brewer.pal(11,'Spectral')[-c(5:7)]
names(rep.col) <- 1:8

## -- treatment Color
treat.col <- c('palegreen3', 'palegreen4', 'skyblue2', 'steelblue','lightsalmon', 'darkorange2') ; names(treat.col) <- names(df.list)[-7]
#plot(1:6, cex=4, col= treat.col, pch=19)
saveRDS(list(asv=colorlist, replicate=rep.col, treat=treat.col), 'Table/03_color_palette_v2.rds')
#######################################################################################
df <- readRDS('Table/03_01_NMDS_all_Sample.rds')

g <- ggplot(df)+
    #geom_path(aes(x=MDS1, y=MDS2, color=treat1))+
    geom_point(aes(x=MDS1, y=MDS2, color=treat1, size=as.numeric(time)))+
    scale_color_manual(values=treat.col)+
    scale_size_area( max_size=3, breaks=seq(1,110, 10))+
    theme_bw()+
    coord_fixed(ratio= max(df$MDS1)/max(df$MDS2))
pdf(sprintf('%s/NMDS.pdf', dir$figdir)); plot(g); dev.off()

a <- prcomp(na.omit(ts[rownames(smllist),]/rowSums(ts[rownames(smllist),])))
df <- cbind(smllist, a$x[rownames(smllist),1:3])

library(animation)
saveGIF(
    for(i in 1:360){
        plot3D::scatter3D(x=df$PC1, y=df$PC2, z=df$PC3, bty = "g",colvar = as.integer(as.factor(df$time)), col = plot3D::gg.col(110),
                          pch=19, phi=0, theta=i)
    },
    movie.name=sprintf("/Users/hiroakifujita/Desktop/Microbiome_TimeSeries/MTS_6/%s/PCA_time.gif", dir$figdir),
    interval=0.01,ani.height=1000,ani.width=1000)

saveGIF(
    for(i in 1:360){
        plot3D::scatter3D(x=df$PC1, y=df$PC2, z=df$PC3, bty = "g",colvar = as.integer(as.factor(df$time)), col = treat.col[order(names(treat.col))],
                                     pch=19, phi=0, theta=i)
    },
    movie.name=sprintf("/Users/hiroakifujita/Desktop/Microbiome_TimeSeries/MTS_6/%s/PCA_media.gif", dir$figdir),
    interval=0.01,ani.height=1000,ani.width=1000)
############################################################################
