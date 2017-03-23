## edgeR figures plotting
## Clear workspace
rm(list=ls(all=T))

## Load libraries
library(ggplot2)
library(gplots)
library(dplyr)
library(devtools)

## Load data

heartVcortex <- read.table("DE_unstranded_heartVScortex.tsv", header = T, sep="\t")
cortexVliver <- read.table("DE_unstranded_cortexVSliver.tsv", header = T, sep="\t")
cortexVlung  <- read.table("DE_unstranded_cortexVSlung.tsv", header = T, sep="\t")

## Select data by strand and omit antisense transcript 8 from analysis
anti.hehi <- subset(heartVcortex, Strand=="+", select=c(GeneID,logFC,PValue,FDR))
anti.hili <- subset(cortexVliver, Strand=="+", select=c(GeneID,logFC,PValue,FDR))
anti.hilu <- subset(cortexVlung, Strand=="+", select=c(GeneID,logFC,PValue,FDR))
anti.hehi$logFC <- anti.hehi$logFC * -1

## Open device for figure output.
pdf(file="figures_output.hg19.pdf", width=12, height=7)

## Set margins for graph
par(mar=c(8.2,4.5,3.8,1))

## par(mfrow=c(1,3))
barplot(anti.hehi$logFC, main = "DE: Cortex VS Heart", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=anti.hehi$GeneID, las = 2, col=1)
abline(h=-1,lty=3)
barplot(anti.hili$logFC, main = "DE: Cortex VS Liver", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=anti.hili$GeneID, las = 2, col=3)
abline(h=-1,lty=3)
barplot(anti.hilu$logFC, main = "DE: Cortex VS Lung", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=anti.hilu$GeneID, las = 2, col=4)
abline(h=-1,lty=3)

## ## par(mfrow=c(1,2))
## barplot(anti.heli$logFC, main = "DE: Heart VS Liver", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=anti.heli$GeneID, las = 2, col=2)
## abline(h=-1,lty=3)
## barplot(anti.helu$logFC, main = "DE: Heart VS Lung", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=anti.helu$GeneID, las = 2, col=5)
## abline(h=-1,lty=3)

## ## par(mfrow=c(1,3))
## barplot(sense.hehi$logFC, main = "DE: Heart VS Cortex", pch = 19, ylab = "logFC", ylim=c(-3,3), names.arg=sense.hehi$GeneID, las = 2, col=1)
## abline(h=3,lty=3)
## barplot(sense.hili$logFC, main = "DE: Cortex VS Liver", pch = 19, ylab = "logFC", ylim=c(-3,3), names.arg=sense.hili$GeneID, las = 2, col=3)
## abline(h=-3,lty=3)
## barplot(sense.hilu$logFC, main = "DE: Cortex VS Lung", pch = 19, ylab = "logFC", ylim=c(-3,3), names.arg=sense.hilu$GeneID, las = 2, col=4)
## abline(h=-3,lty=3)

## ## par(mfrow=c(1,2))
## barplot(sense.heli$logFC, main = "DE: Heart VS Liver", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=sense.heli$GeneID, las = 2, col=2)
## abline(h=-1,lty=3)
## barplot(sense.helu$logFC, main = "DE: Heart VS Lung", pch = 19, ylab = "logFC", ylim=c(-2,0), names.arg=sense.helu$GeneID, las = 2, col=5)
## abline(h=-1,lty=3)

dev.off()
