##Plot unstranded human data comparing tissues (cortex to heart, liver, and lung)
##Clear workspace
rm(list=ls(all=T))

##Load libraries
library(dplyr)
library(ggplot2)
library(devtools)

pdf(file="tissue.comparison.hg19.pdf", useDingbats=FALSE)

exon.unstranded <- read.table("unstranded_exon.data.tsv", sep="\t", header=TRUE)

fc.exon <- ggplot(data=exon.unstranded, aes(x=Sample, y=logFC, fill=GeneID)) +
    geom_bar(colour="black", width=0.8, stat="identity") +
    geom_hline(aes(yintercept=0))

fc.exon + scale_x_discrete(name="") +
    scale_y_continuous(name="Expression: log2(FC)",limits=c(-6,1)) + 
    scale_fill_manual(values="grey25") +
    annotate("text",x=1,y=-5.40,label="***,0.00", size=5) +
    annotate("text",x=2,y=-5.14,label="***,0.00", size=5) +
    annotate("text",x=3,y=-5.83,label="***,0.00", size=5) +
    theme(legend.title=element_blank(),
          legend.position=c(0.88,0.95),
          legend.text=element_text(size=16),
          axis.text.x=element_text(face="bold", colour="black", size=16),
          axis.title.y=element_text(face="bold", size=18),
          axis.text.y=element_text(size=16,colour="black"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill="white"))

dev.off()
