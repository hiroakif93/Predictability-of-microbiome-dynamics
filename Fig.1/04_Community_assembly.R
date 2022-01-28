############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot', 'extrafont','RColorBrewer', 'scales', 'ggsnippets'))

# -- Create directory to save
dir <- make.dir('Fig.1/Community assembly')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/Abruptness.rds')

taxa <- as.matrix(readRDS('Table/Taxa_list.rds'))
bin_info <- unique(read.table('Table/bin_summary.txt', sep='\t', header=T)) 

color <- readRDS('Table/color_palette.rds')

## -- ASV X_0003 showed same dynamics with MAG
taxa['X_0003', 2:8] <- as.character(as.matrix(bin_info)[22,4:10])
############################################################################

glist1 <- c()
names(dlist) <- gsub('medium ', 'Medium-', names(dlist))
sml$treat1 <- gsub('medium ', 'Medium-', sml$treat1)

for(i in names(dlist)[-7]){ #i=names(dlist)[2]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]
    rel.ts <- ts/rowSums(ts)
    
    ## ======================================= ##
    ## -- Visualizing community assembly
    l='Genus'
        agg.ts <- t(Taxa.mat(ts[rownames(smlsub),], taxa, l))
        
        lf <- gather(cbind(smlsub, agg.ts), 
                     key, value, -c(1:(ncol(smlsub))))
        
        lf$key <- factor(lf$key, levels=names(sort(colSums(agg.ts))) )
        
        col <- color[[1]][[l]][colnames(agg.ts)]
        if(l=='Genus' & any(is.na(col))) {names(col)[is.na(col)] <- taxa['X_0003', 7] ;col[is.na(col)] <- 'yellow'}
        if(any(names(col)%in%taxa['X_0003', c(1,5:7)])) col[names(col)%in%taxa['X_0003', c(1,5:7)]] <- 'yellow'

        g1 <- ggplot(lf)+
            geom_area(aes(x=as.numeric(time), y=value, fill=key), color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
            facet_wrap(~replicate.id,scales='free_y',ncol=1, strip.position='left')+
            scale_fill_manual( values=col ) + 
            theme_bw(base_size=15)+
            labs(x= "Day", subtitle=i, y= "Number of DNA coies")+
            guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))+
            theme(text=element_text(family='Arial'),
                  strip.placement = 'outside')+
            scale_x_continuous(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0), position = "right")
        glist1[[i]] <- g1
       
}    

ggsave(plot=plot_grid(plotlist=glist1, nrow=2, byrow=FALSE),
       filename=sprintf('%s/Absolute abundance of community structure (Ex.Fig2).pdf',
       dir$figdir, gsub('/', '_', i)),h=21, w=15)
############################################################################
glist1 <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[4]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]
    rel.ts <- ts/rowSums(ts)
    a <- smlsub[which(smlsub[,'abruptness_rel_tw5_tp1']>0.5), ]

    ## ======================================= ##
    ## -- Visualizing community assembly
    
    agg.ts <- t(Taxa.mat(ts[rownames(smlsub),], taxa, l))
    
    lf <- gather(cbind(smlsub, agg.ts), 
                 key, value, -c(1:(ncol(smlsub))))
    
    lf$key <- factor(lf$key, levels=names(sort(colSums(agg.ts))) )
    
    col <- color[[1]][[l]][colnames(agg.ts)]
    if(l=='Genus' & any(is.na(col))) {names(col)[is.na(col)] <- taxa['X_0003', 7] ;col[is.na(col)] <- 'yellow'}
    if(any(names(col)%in%taxa['X_0003', c(1,5:7)])) col[names(col)%in%taxa['X_0003', c(1,5:7)]] <- 'yellow'
    
    g1 <- ggplot(lf)+
        geom_area(aes(x=as.numeric(time), y=value, fill=key), position='fill', color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
        geom_line(aes(x=as.numeric(time), y=abruptness_rel_tw5_tp1), size=1.5)+
        geom_line(aes(x=as.numeric(time), y=abruptness_rel_tw5_tp1),
                  size=1.3, color='olivedrab1')+
        geom_hline(yintercept=0.5, linetype=2)+
        facet_wrap(~replicate.id,scales='free_y',ncol=1, strip.position='left')+
        scale_fill_manual( values=col ) + 
        theme_bw(base_size=15)+
        labs(x= "Day", subtitle=i, y= "Relative abundance")+
        guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))+
        theme(text=element_text(family='Arial'),
              strip.placement = 'outside')+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0), position = "right")
    glist1[[i]] <- g1
    
}    

ggsave(plot=plot_grid(plotlist=glist1, nrow=2, byrow=FALSE),
       filename=sprintf('%s/Relative abundance of community structure (Ex.Fig3).pdf', dir$figdir, gsub('/', '_', i)),h=21, w=15)
       
#######################################################################################
     
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
	saveRDS(df, 'Table/NMDS_all_Sample.rds')
	saveRDS(nmds, 'Table/NMDS_stat.rds')
}else{
	df <- readRDS('Table/NMDS_all_Sample.rds')
}

g <- ggplot(df)+
    #geom_path(aes(x=MDS1, y=MDS2, color=treat1))+
    geom_point(aes(x=MDS1, y=MDS2, color=treat1, fill=after_scale(alpha(color, 0.7)), size=as.numeric(time)),
               alpha=0.7, shape=21)+
    scale_color_manual(values= color[[3]], guide=guide_legend(ncol=2, direction='horizontal'))+
    scale_size_area( max_size=5, breaks=c(1, 50, 110), 
                     guide=guide_legend(ncol=1, title.position = 'right'))+
    theme_bw(base_size=35)+
    theme(text=element_text(family='Arial'),
          legend.position='bottom',
          plot.margin= unit(c(1, 3, 1, 1), "lines"))+
    labs(x='Axis 1', y='Axis 2', size='', color='')+
    coord_fixed(ratio= (max(df$MDS1)-min(df$MDS1))/(max(df$MDS2)-min(df$MDS2)))

ggsave(plot=g, file=sprintf('%s/NMDS (Fig.1-d).pdf', dir$figdir),
       w=10, h=15)
