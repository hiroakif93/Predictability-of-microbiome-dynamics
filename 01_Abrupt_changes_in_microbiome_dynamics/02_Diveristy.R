############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('../')
#### 
############################################################################

## -- Loading Function and Library
source('functions/functions.R')
load.lib( c('vegan', 'ggplot2', 'tidyr', 'cowplot', 'RColorBrewer', 'extrafont','doParallel', 'scales'))

# -- Create directory to save
dir <- make.dir('01_Abrupt_changes_in_microbiome_dynamics/Diversity')

# -- Load data table
dlist  <- readRDS('Table/matrixList.rds')
sml <- readRDS('Table/sample_info.rds')
saveRDS(sml, 'Table/sample_info.rds')

taxa <- readRDS('Table/Taxa_list.rds')
color <- readRDS('Table/color_palette.rds')

############################################################################

BoxMeanQuant <- function(x) {
    v <- c(quantile(x, 0.25), quantile(x, 0.25), median(x), quantile(x, 0.75), quantile(x, 0.75))
    names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
    v
}


smlsum <- NULL
for(i in names(dlist)[-7]){ 
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
	
    ts.sub <- na.omit(dlist[[i]]) 
    smlsub <- sml[rownames(ts.sub),]
    
    ## ========================================================== ##
    ## -- Alpha diveristy (Shannon's H')
    smlsub$shannon <- apply(ts.sub, 1, diversity, index='shannon')
	
	## ========================================================== ##
 	## -- Beta diversity against replicate
 	
 	smlsub$bray_against_replicate <- NA; smlsub$beta.q0_against_replicate <- NA;
	for(n in 1:110){#n=50
		
		smlday <- smlsub[ which(as.numeric(smlsub$time)==n), ]
		
		bray <- colMeans( as.matrix(vegdist(ts.sub[rownames(smlday), ], method='bray')) )

		smlsub[names(bray), ]$bray_against_replicate <- bray
	}
	
	## ========================================================== ##
	
	smlsum <- rbind(smlsum, smlsub)
}

saveRDS(smlsum[,9:10], 'Table/01_Diversity.rds')
sml <- smlsum
############################################################################
## -- Mean coefficient of variation of absolute abundance
names(dlist) <- gsub('medium ', 'Medium-', names(dlist))
sml$treat1 <- gsub('medium ', 'Medium-', sml$treat1)

sml$treat1 <- factor(sml$treat1, levels=names(dlist)[-7])

cvs <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[1]
    
    ts <- dlist[[ i ]]
    smlsub <- sml[sml$treat1==i, ]
    
    replicate <- c()
    for(r in 1:8){
        
        smlrep <- smlsub[smlsub$replicate.id==r, ]
        tssub <- ts[rownames(smlrep), ]

        cv <- apply(tssub[, colSums(tssub>0)>10], 2, function(x) sd(x)/mean(x) )
        
        replicate <- rbind(replicate, cbind(replicate=r, cv))
    }
    
    cvs  <- rbind(cvs, cbind(treat=i, as.data.frame(replicate)))

}

cvg <- ggplot(cvs, aes(x=treat, y=cv))+
        geom_jitter(aes(x=treat, y=cv, color=as.factor(replicate)))+
        stat_summary(fun.data = BoxMeanQuant, geom = "boxplot",width=0.7,show.legend=F, size=0.3, alpha=0.8 )+
        scale_color_manual(values=color$replicate)+
        theme_bw(base_size=15) +
        labs(x='', y="Mean coefficient of variation\nof absolute abundance", color='Replicate\ncommunity')+
        theme(axis.text.x=element_text(angle=45, vjust=0.98, hjust=1),
              text=element_text(family='Arial'),
              panel.grid=element_blank(),
              legend.title=element_text(size=10),
              legend.text=element_text(size=7),
              legend.key.size = unit(0.7,"line"))

ggsave(plot= cvg,
       file=sprintf('%s/Mean coefficient of variation.tiff', dir$figdir),
       w=6.5, h=6)

############################################################################

alpha <- ggplot(sml, aes(x=treat1, y=shannon))+
    geom_jitter(aes(x=treat1, y=shannon, color=as.factor(replicate.id)))+
    stat_summary(fun.data = BoxMeanQuant, geom = "boxplot",width=0.7,show.legend=F, size=0.3, alpha=0.8 )+
    scale_color_manual(values=color$replicate)+
    theme_bw(base_size=15) +
    labs(y="Alpha diveristy\n(Shannon's H')", x="", color='Replicate\ncommunity')+
    theme(axis.text.x=element_text(angle=45, vjust=0.98, hjust=1),
          text=element_text(family='Arial'),
          panel.grid=element_blank(),
          legend.key.size = unit(0.1,"line"))

ggsave(plot=alpha,
       file=sprintf('%s/Alpha diveristy.tiff', dir$figdir),
       w=6.5, h=6)

# -- Visualize alpha diversity dynamics (species richness)
g <- ggplot(data= sml, aes(x=as.numeric(time), y= shannon,group = replicate.id,color=as.character(replicate.id))) +
    geom_point(alpha=0.8,na.rm=T) + 
    #geom_smooth(method=lm,formula = y~splines::bs(x,3),alpha=0.4, size=1,se=F) + 
    scale_colour_manual(values = brewer.pal(11,'Spectral')[-c(5:7)])+
    facet_wrap(~treat1, dir='h', ncol=2) + 
    labs(x='Day',y="Alpha diveristy\n(Shannon's H')",color='Replicate\ncommunity') + 
    theme_bw(base_size= 20, base_family='') +
    theme(text=element_text(family='Arial'),
          strip.background = element_blank(),
          panel.grid =  element_line(size=0.5, colour = "azure2"),
          panel.grid.minor = element_line(size=0.5, colour = "azure2"),
          panel.border = element_rect(fill=NA,colour = "black"),
          #plot.title = element_text(size = 40,  face = "bold"),
          legend.position='bottom')

ggsave(plot=g,
       file=sprintf('%s/Time series of Alpha diveristy.tiff', dir$figdir),
       width=10,heigh=13)
       
          
############################################################################
# -- Beta diversity of each treatment

statsum <- gglist <- c()
for( i in names(dlist)[-7]){ # for test : i=names(dlist)[1]
    
    ## ========================================================== ##
    ## -- Extracting one treatment matrix
    
    ts <- na.omit(dlist[[i]]) 
    smlsub <- sml[which(sml$treat1==i),]
    smlsub <- smlsub[order(smlsub$time), ]
    
    ## ========================================================== ##
    ## -- Calculate Bray-Curtis dissimilarity
    
    lf.asv <- lf.fun <- beta<- cor.r <- c() ; fn.bg <- ab.bg <- list()
    for(r in 1:8){ # r=1
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        # -- Data compile
        smlrep <- smlsub[smlsub$replicate.id == r, ]
        tssub <- ts[rownames(smlrep), ] ; 
        
        ## -- Absolute value / relative value
        mat.list <- list(tssub, tssub/rowSums(tssub))
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        # -- Calculate Bray-Curtis dissimilarity
        sample <- rownames(smlrep)
        beta.list <- lapply(mat.list, function(x){ # for test ; x=mat.list[[1]]
            beta <- as.matrix( vegdist( x, method='bray') )
            beta.mat <- t(as.data.frame(beta)[sample,] )
            beta.mat <- as.data.frame(beta.mat)[sample,] 
            colnames(beta.mat) <- smlrep[,'time']
            beta.mat2 <- as.data.frame(cbind(x=smlrep$time,replicate=paste('Replicate', r), beta.mat)) 
            
        })
        
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        # -- Mean beta diversity
        beta.mean <- do.call(cbind,lapply(beta.list, function(x){colMeans(x[,-c(1:2)], na.rm=TRUE)}))
        colnames(beta.mean) <- c('abs.betaMean.in.replicate','rel.betaMean.in.replicate')

        # -- Compile to long format for ggplot
        lf <- lapply(beta.list, function(l){ gather(l, key=key, value=value, -c(1:2)) } )  	
        lf.asv <- rbind(lf.asv, lf[[1]]) ; lf.fun <- rbind(lf.fun, lf[[2]]) 
        
        # -- Calculate Bray-Curtis dissimirarity witin froom one to 'beta.lag' day
        betalag.tp1 <- lapply(beta.list, function(x){# for test ; x= beta.list[[1]]
            a <- x[,-c(1,2)]
            c(NA, diag(as.matrix(a[1:(nrow(a)-1), 2:ncol(a)]))) })
        betalag.asv <- betalag.tp1[[1]];   betalag.fun <- betalag.tp1[[2]]  		
        beta <- rbind(beta, cbind(smlrep, asv.beta_lag =betalag.asv, function.beta_lag =betalag.fun, beta.mean))
        
    }
    
    statsum <- rbind(statsum, beta)
    
    
    ## -------------------------------------------------------------------------------------- ##
    
    ## -- Visualize results
    lf.asv$x <- factor(lf.asv$x, levels=unique(lf.asv$x)) ; lf.fun$x <- factor(lf.fun$x, levels=unique(lf.fun$x))
    gg.ab <- ggplot() +
        geom_tile(data=lf.asv, aes(x=as.numeric( gsub('Day_','', x)), y=as.numeric( gsub('Day_','', key)),fill=as.numeric(value)))+
        scale_fill_gradientn(colors=brewer.pal(8,'RdGy'),na.value='black') +
        labs(y='Day', x='Day',fill='score')+
        scale_y_continuous(expand=c(0,0), breaks=c(1,110)) + scale_x_continuous(expand=c(0,0), breaks=c(110)) +
        facet_wrap(~replicate) + 
        theme_minimal(base_size=13)  +
        theme(panel.grid = element_blank(),
              legend.position='bottom',
              legend.justification=c(1,0),
              legend.title=element_text(vjust=0.8),
              legend.key.width=unit(1,'cm'),
              legend.key.heigh=unit(0.3,'cm'),
              axis.ticks = element_blank(),
              text=element_text(family='Arial'))+
        coord_equal()  
    
    gg.fn <- ggplot() +
        geom_tile(data= lf.fun, aes(x=as.numeric( gsub('Day_','', x)), y=as.numeric( gsub('Day_','', key)),fill=as.numeric(value) ))+
        scale_fill_gradientn(colors=brewer.pal(8,'RdGy'),na.value='black') +
        labs(y='Day', x='Day',fill='score')+
        scale_y_continuous(expand=c(0,0), breaks=c(1,110)) + scale_x_continuous(expand=c(0,0), breaks=c(110)) +
        facet_wrap(~replicate) + 
        theme_minimal(base_size=13)  +
        theme(panel.grid = element_blank(),
              legend.position='bottom',
              legend.justification=c(1,0),
              legend.title=element_text(vjust=0.8),
              legend.key.width=unit(1,'cm'),
              legend.key.heigh =unit(0.3,'cm'),
              axis.ticks = element_blank(),
              text=element_text(family='Arial'))+
        coord_equal() 
    gg.fn2 <- ggplot() +
        geom_tile(data= lf.fun, aes(x=as.numeric( gsub('Day_','', x)), y=as.numeric( gsub('Day_','', key)),fill=as.numeric(value) ), show.legend=F)+
        scale_fill_gradientn(colors=brewer.pal(8,'RdGy'),na.value='black') +
        labs(x='Day', 
             y='Day',fill='score', subtitle=i)+
        scale_y_continuous(expand=c(0,0), breaks=c(1,110)) + scale_x_continuous(expand=c(0,0), breaks=c(110)) +
        facet_wrap(~replicate) + 
        theme_minimal(base_size=13)  +
        theme(panel.grid = element_blank(),
              legend.position='bottom',
              legend.justification=c(1,0),
              legend.title=element_text(vjust=0.8),
              legend.key.width=unit(1,'cm'),
              legend.key.heigh =unit(0.3,'cm'),
              axis.ticks = element_blank(),
              text=element_text(family='Arial'))+
        coord_equal() 
    gglist[[i]] <- gg.fn2


}

ggsave(plot=plot_grid(plotlist=gglist, byrow=FALSE),
       filename=sprintf('%s/Bray-Curtis dissimilarity of community structure between time points.tiff', dir$figdir),
       width=10,heigh=8)

## -------------------------------------------------------------------------------------- ##
# -- Beta diversity dynamics among replicate
betadiv <- c()
for(i in names(dlist)[-7]){ #i=names(dlist)[1]
	
	## ========================================================== ##
    ## -- Extracting one treatment matrix
    ts <- dlist[[i]]
    smlsub <- sml[sml$treat1==i, ] 
    tarel <- ts/rowSums(ts)

	tmp <- (ts/rowSums(ts))[ rownames(smlsub[smlsub$replicate.id==5,]) ,]
	
    ## ========================================================== ##
	## -- indices
	
	## -- Against replicate
	beta <- matrix(NA, ncol=1, nrow=110,
				   dimnames=list(unique(smlsub$time), 'bray_asv'))
	for(d in as.character(unique(smlsub$time))){ #d=as.character(unique(smlsub$time))[7]
		
		smlday <- smlsub[smlsub$time==d, ]
		beta[d,1] <- mean( vegdist(na.omit(tarel[rownames(smlday), ]), method='bray', na.rm=TRUE), na.rm=TRUE)

	}

	## ========================================================== ##
	## -- visualizing
	
	df <- cbind(treat=i, time=1:110, as.data.frame(beta) )
	betadiv <- rbind(betadiv, df)
}

betadiv$treat <- factor(betadiv$treat, sort(as.character(unique(betadiv$treat))) )
betadiv2 <- cbind(betadiv, group=do.call(rbind, strsplit(as.character(betadiv$treat), '/')))
betadiv2 $group.2 <- gsub('me', 'Me', betadiv2 $group.2)

g1 <- ggplot(betadiv2)+
	  geom_line(aes(x=time, y=bray_asv, group=group.1, color =group.1), size=2)+
	  facet_wrap(~ group.2)+
	  theme_bw(base_size=26)+
	  scale_color_manual(values=c('brown', 'skyblue3'), guide=guide_legend(ncol=2, title.position='left'))+
	  labs(color='Inoculum', x='Day', y='Bray-Curtis dissimilarity\namong replicate')+
	  theme(text=element_text(family='Arial'),
          panel.grid =  element_line(size=0.5, colour = "azure2"),
          panel.grid.minor = element_line(size=0.5, colour = "azure2"),
          panel.border = element_rect(fill=NA,colour = "black"),
          legend.position=c(1, 0),
          legend.justification = c(1,2),
          strip.background = element_blank(),
          plot.margin= unit(c(1, 1, 7, 1), "lines"))

ggsave(plot=g1, file=sprintf('%s/Bray-Curtis dissimilarity among replicate.tiff', dir$figdir), h=8, w=21) 