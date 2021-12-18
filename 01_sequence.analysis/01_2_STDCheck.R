############################################################################
####
#### R script for Fujita (2019)
####
#### Check standard DNA quolity
#### 2019.10.31 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('..')
#### setwd('~/Desktop/Microbiome_TimeSeries/MTS_5')
############################################################################

# -- Create directory to save
wd.dir <- '01_sequence.analysis'
save.dir <- sprintf('%s/01_output', wd.dir)

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load workspace and functions
list <- readRDS('Table/for_STD_check.rds')
taxa.print <- list$asv.table
seqtab <- list$seqtab

source(sprintf('%s/function_refrence/F1_HelperFunctions.R',wd.dir))

# -- Load library and functions
library(ggplot2); packageVersion("ggplot2") 	# v3.2.1
library(cowplot); packageVersion("cowplot") 	# v1.0.0
library(reshape2); packageVersion("reshape2")	# v1.4.3
library(stringr); packageVersion("stringr") 	# v1.4.0

# -- Nnumber of copies of standard DNA
std.copy.n <- c( 0.1, 0.05, 0.02, 0.01,0.005)*(6.02*10^13/1000000)/1000

############################################################################

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
new.seqtab <- as.data.frame(seqtab[,-n.std.seq]) # Make seq table without standard DNA
new.seqtab2 <- new.seqtab

colnames(new.seqtab2) <- taxa.print[-n.std.seq, "Phylum"]

# --  Visualize regression results
# --  1. R2 value distribution
g1 <- ggplot(data.frame(values = r2.summary), aes(x = values))
g1 <- g1 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2) #+ xlim(0.9,1) + ylim(0,13)
g1 <- g1 + xlab(expression(paste("R"^{2}, " values"))) + ylab("Count") #+ scale_x_continuous(breaks=seq(0.8,1,by=0.05))
g1 <- PlotStyle(g1)

# --  2. Slope distribution
g2 <- ggplot(data.frame(values = coef.summary), aes(x = values))
g2 <- g2 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2)
g2 <- g2 + xlab("Regression slope") + ylab("Count")
g2 <- PlotStyle(g2)

# --  3. Regression examples
max.slope <- as.numeric(c(new.std.table[which.max(coef.summary),], 0))
med.slope <- as.numeric(c(new.std.table[which.min(abs(coef.summary - median(coef.summary))),], 0))
min.slope <- as.numeric(c(new.std.table[which.min(coef.summary),], 0))
slope.summary <- melt(data.frame(copy =c(std.copy.n, 0),
                                 max_slope = max.slope,
                                 med_slope = med.slope,
                                 min_slope = min.slope), id.vars = "copy")
g3 <- ggplot(slope.summary, aes(x = copy, y = value, group = variable, colour = variable))
g3 <- g3 + geom_point(size = 2) + scale_color_manual(name = "Regression slope", values = c("red3", "darkred", "black"))
g3 <- g3 + geom_smooth(method = "lm", size = 0.5, se = F)+theme_bw(base_size=26)
	  
g3 <- g3 + xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads")
g3 <- PlotStyle(g3) + theme(text=element_text(family='Arial'),legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))

ggsave(plot= g3, w=5, h=3,
		filename=sprintf("%s/STD_slope_exfig.tiff", save.dir))

# --  4. Read visualization
new.seqtab2$sample <- rownames(new.seqtab2)

new.seqtab2[]
seqtab.plot <- melt(new.seqtab2, id.vars = "sample")
seqtab.plot <- seqtab.plot[seqtab.plot$value > 0,]
seqtab.plot$value2 <- seqtab.plot$value + 0.5
std.table.plot <- data.frame(sample = rownames(new.seqtab2),
                             max = apply(new.std.table, 1, max),
                             min = apply(new.std.table, 1, min) + 0.5) # data.frame for overwrite standard DNA reads

# --  Box plot + jitter plot version
g4 <- ggplot(seqtab.plot, aes(x = sample, y = value2))
g4 <- g4 + geom_boxplot(shape = 16, alpha = 0.5)
g4 <- g4 + geom_jitter(shape = 16, size = 1.5, alpha = 0.5, position = position_jitter(0.2))
g4 <- g4 + scale_y_log10() + ylab("Sequence reads") + xlab("Sample ID")
g4 <- PlotStyle(g4) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Overwrite standard DNA reads
g4 <- g4 + geom_segment(data = std.table.plot, aes(x = sample, xend = sample, y = min, yend = max),
                        colour = "red3", size = 1.5, alpha = 0.5)
g4 <- g4 + ggtitle("QMiSeq sequence analysis demo")

# -- Summarize visualization
top.row <- plot_grid(g1, g2, g3, ncol = 3, align = "h", labels = c("A", "B", "C"), rel_widths = c(1,1,1.3))
Fig.std <- plot_grid(top.row, g4, ncol = 1, labels = c("", "D"))

############################################################################

# -- Save figure
pdf(sprintf("%s/STD_Check.pdf", save.dir),heigh=10,width=15)
Fig.std
dev.off()

# Save and output results
save.image(sprintf("%s/01_2_output.RData",save.dir))

############################################################################
