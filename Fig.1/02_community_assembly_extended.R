############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot', 'extrafont','RColorBrewer', 'scales', 'ggsnippets'))

# -- Create directory to save
dir <- make.dir('99_Figure_for_MS/Fig1')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS('Table/03_matrixList.rds')
sml <- readRDS('Table/03_06_sample_info.rds')

taxa <- as.matrix(readRDS('Table/Taxa_list_02.rds'))
bin_info <- unique(read.table('Table/08_metagenome/bin_summary.txt', sep='\t', header=T)) ; dim(bin_info)

color <- readRDS('Table/99_color_palette_v2.rds')

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
        
        col <- color[[l]][colnames(agg.ts)]
        if(l=='Genus' & any(is.na(col))) {names(col)[is.na(col)] <- taxa['X_0003', 7] ;col[is.na(col)] <- 'yellow'}
        if(any(names(col)%in%taxa['X_0003', c(1,5:7)])) col[names(col)%in%taxa['X_0003', c(1,5:7)]] <- 'yellow'

        g1 <- ggplot(lf)+
            geom_area(aes(x=as.numeric(time), y=value, fill=key), color='grey30', stat='identity', show.legend=FALSE, size=0.1)+
            geom_vline(xintercept=83)+
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
pdf(sprintf('%s/community_assembly.pdf', dir$figdir, gsub('/', '_', i)),h=21, w=15)
plot(plot_grid(plotlist=glist1, nrow=2, byrow=FALSE))
dev.off()

ggsave(plot=plot_grid(plotlist=glist1, nrow=2, byrow=FALSE),
       filename=sprintf('%s/community_assembly.tiff', dir$figdir, gsub('/', '_', i)),h=21, w=15)
############################################################################
glist1 <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[4]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]
    rel.ts <- ts/rowSums(ts)
    a <- smlsub[which(smlsub[,'abruptness_rel_tw5_tp1']>0.5), ]
    
    a[,1:10]
    ## ======================================= ##
    ## -- Visualizing community assembly
    
    agg.ts <- t(Taxa.mat(ts[rownames(smlsub),], taxa, l))
    
    lf <- gather(cbind(smlsub, agg.ts), 
                 key, value, -c(1:(ncol(smlsub))))
    
    lf$key <- factor(lf$key, levels=names(sort(colSums(agg.ts))) )
    
    col <- color[[l]][colnames(agg.ts)]
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
        labs(x= "Day", subtitle=i, y= "Number of DNA coies")+
        guides(fill=guide_legend(title='',ncol=3, reverse=TRUE))+
        theme(text=element_text(family='Arial'),
              strip.placement = 'outside')+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0), position = "right")
    glist1[[i]] <- g1
    
}    
pdf(sprintf('%s/community_assembly_relative.pdf', dir$figdir, gsub('/', '_', i)),h=21, w=15)
plot(plot_grid(plotlist=glist1, nrow=2, byrow=FALSE))
dev.off()


ggsave(plot=plot_grid(plotlist=glist1, nrow=2, byrow=FALSE),
       filename=sprintf('%s/community_assembly_relative.tiff', dir$figdir, gsub('/', '_', i)),h=21, w=15)