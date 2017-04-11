#Ballgown Analysis for Thesis
#Clear workspace
rm(list=ls(all=T))

#Load libraries
library(ballgown)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(genefilter)
library(dplyr)
library(devtools)

#Load phenotype data to generate ballgown object
pheno_unstranded <- read.csv("unstranded.hg19.csv")
pheno_stranded   <- read.csv("stranded.hg19.csv")
pheno_neurons    <- read.csv("neurons.hg19.csv")

## ##Generated ballgown object and save
## bg_unstranded <- ballgown(samples=as.vector(pheno_unstranded$path),pData=pheno_unstranded)
## bg_stranded   <- ballgown(samples=as.vector(pheno_stranded$path),pData=pheno_stranded)
## bg_neurons    <- ballgown(samples=as.vector(pheno_neurons$path),pData=pheno_neurons)

## save(bg_neurons, file='./bg_neurons.rda')
## save(bg_stranded, file='./bg_stranded.rda')
## save(bg_unstranded, file='./bg_unstranded.rda')

#Load ballgown variable for script
load("bg_neurons.rda")   #40931 transcripts and 10 samples
load("bg_stranded.rda")  #40931 transcripts and 6 samples
load("bg_unstranded.rda")#40931 transcripts and 24 samples

#Filter data
bg_unstranded.filt <- subset(bg_unstranded, "rowVars(texpr(bg_unstranded)) > 1", genomesubset=TRUE)
bg_stranded.filt   <- subset(bg_stranded, "rowVars(texpr(bg_stranded)) > 1", genomesubset=TRUE)
bg_neurons.filt    <- subset(bg_neurons, "rowVars(texpr(bg_neurons)) > 1", genomesubset=TRUE)

##For neurons only (two sample comparison)
#Statistical significant test, add transcript names, arrange by p.value, and write to .csv file
## bg_neurons.table    <- texpr(bg_neurons, 'all')
## bg_neurons.names    <- unique(bg_neurons.table[,c(1,6)])
## neurons_transcripts <- stattest(bg_neurons, feature="transcript", covariate="type", getFC=T, meas="FPKM")
## neurons_transcripts <- merge(neurons_transcripts, bg_neurons.names, by.x=c("id"),by.y=c("t_id"))
## neurons_transcripts <- arrange(neurons_transcripts,pval)
## write.csv(neurons_transcripts, "neurons.transcript_results.csv", row.names=F)

bg_unstranded.table.filt <- texpr(bg_unstranded.filt, 'all')
bg_unstranded.table      <- texpr(bg_unstranded, 'all')
bg_unstranded.names      <- unique(bg_unstranded.table.filt[,c(1,6)])
bg_unstranded.names.all  <- unique(bg_unstranded.table[,c(1,6)])
unstranded_transcripts   <- stattest(bg_unstranded.filt, feature="transcript", covariate="type", getFC=F, meas="FPKM")
unstranded_transcripts.all <- stattest(bg_unstranded, feature="transcript", covariate="type", getFC=F, meas="FPKM")
unstranded_transcripts   <- merge(unstranded_transcripts, bg_unstranded.names, by.x=c("id"),by.y=c("t_id"))
unstranded_transcripts.all <- merge(unstranded_transcripts.all, bg_unstranded.names.all, by.x=c("id"),by.y=c("t_id"))
unstranded_transcripts   <- arrange(unstranded_transcripts,pval)
## write.csv(unstranded_transcripts, "unstranded.transcript_results.filtered.csv", row.names=F)
## write.csv(unstranded_transcripts.all, "unstranded.transcript_results.csv", row.names=F)

## bg_stranded.table.filt <- texpr(bg_stranded.filt, 'all')
## bg_stranded.names      <- unique(bg_stranded.table.filt[,c(1,6)])
## stranded_transcripts   <- stattest(bg_stranded.filt, feature="transcript", covariate="type", getFC=F, meas="FPKM")
## stranded_transcripts   <- merge(stranded_transcripts, bg_stranded.names, by.x=c("id"),by.y=c("t_id"))
## stranded_transcripts   <- arrange(stranded_transcripts,pval)
## write.csv(stranded_transcripts, "stranded.transcript_results.filtered.csv", row.names=F)

#Assign values to possible isoforms for UBE3A and UBE3A antisense
#UBE3A antisense transcripts
iso1.1  <- 15947 #MSTRG.178.3
iso1.2  <- 15944 #MSTRG.178.2
iso1.3  <- 15943 #MSTRG.178.1
iso1.4  <- 15950 #MSTRG.178.5
iso1.5  <- 15949 #MSTRG.178.4

iso2.1  <- 15955 #MSTRG.178.7
iso2.2  <- 15957 #MSTRG.178.8
iso2.3  <- 15966 #MSTRG.178.10

iso3.1  <- 15970 #MSTRG.178.13
iso3.2  <- 15974 #MSTRG.178.14

circRNA <- 16029 #MSTRG.178.26

#UBE3A transcripts	
UBE3A.iso1.1 <- 16021	#MSTRG.222.6
UBE3A.iso1.2 <- 16026	#MSTRG.222.11

UBE3A.iso2.1 <- 16016	#MSTRG.222.1
UBE3A.iso2.2 <- 16019	#MSTRG.222.4
UBE3A.iso2.3 <- 16020	#MSTRG.222.5
UBE3A.iso2.4 <- 16027	#MSTRG.222.13
UBE3A.iso2.5 <- 16028	#MSTRG.222.14

UBE3A.iso4.1 <- 16017	#MSTRG.222.2
UBE3A.iso4.2 <- 16018	#MSTRG.222.3

UBE3A.iso5   <- 16025	#MSTRG.222.12
UBE3A.iso6   <- 16015	#MSTRG.222.7

#Boxplot comparing expression of the different possible isoform assemblies for antisense UBE3A and UBE3A
pdf(file="ballgown_unstranded_output.pdf")

fpkm.unstranded.all <- texpr(bg_unstranded, meas="FPKM")
fpkm.unstranded.all <- log2(fpkm.unstranded.all+1)

boxplot(fpkm.unstranded.all[iso1.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso1.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso1.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso1.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso1.3,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso1.3,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso1.4,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.4]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso1.4,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso1.5,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso1.5,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[iso2.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso2.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso2.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso2.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso2.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso2.3,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso2.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso2.3,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[iso3.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso3.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso3.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[iso3.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[iso3.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[iso3.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[circRNA,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[circRNA]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[circRNA,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

#UBE3A transcripts
boxplot(fpkm.unstranded.all[UBE3A.iso1.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso1.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso1.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso1.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso1.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[UBE3A.iso2.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso2.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso2.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso2.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso2.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso2.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso2.3,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso2.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso2.3,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso2.4,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso2.4]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso2.4,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso2.5,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso2.5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso2.5,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[UBE3A.iso4.1,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso4.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso4.1,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))
boxplot(fpkm.unstranded.all[UBE3A.iso4.2,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso4.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso4.2,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[UBE3A.iso5,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso5,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

boxplot(fpkm.unstranded.all[UBE3A.iso6,] ~ pheno_unstranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_unstranded)[UBE3A.iso6]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.unstranded.all[UBE3A.iso6,] ~ jitter(as.numeric(pheno_unstranded$type)), col=as.numeric(pheno_unstranded$type))

#Average expression of transcripts
plotMeans("MSTRG.222",bg_unstranded, groupvar="type", legend=T)
plotMeans("MSTRG.178",bg_unstranded, groupvar="type", legend=T)

dev.off()

#Boxplot comparing expression of the different possible isoform assemblies for antisense UBE3A and UBE3A
pdf(file="ballgown_stranded_output.pdf")

fpkm.stranded.all <- texpr(bg_stranded, meas="FPKM")
fpkm.stranded.all <- log2(fpkm.stranded.all+1)

boxplot(fpkm.stranded.all[iso1.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso1.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso1.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso1.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso1.3,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.3]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso1.3,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso1.4,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.4]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso1.4,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso1.5,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.5]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso1.5,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[iso2.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso1.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso2.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso2.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso2.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso2.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso2.3,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso2.3]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso2.3,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[iso3.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso3.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso3.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[iso3.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[iso3.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[iso3.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[circRNA,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[circRNA]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[circRNA,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

#UBE3A transcripts
boxplot(fpkm.stranded.all[UBE3A.iso1.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso1.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso1.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso1.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso1.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso1.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[UBE3A.iso2.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso2.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso2.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso2.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso2.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso2.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso2.3,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso2.3]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso2.3,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso2.4,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso2.4]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso2.4,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso2.5,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso2.5]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso2.5,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[UBE3A.iso4.1,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso4.1]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso4.1,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))
boxplot(fpkm.stranded.all[UBE3A.iso4.2,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso4.2]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso4.2,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[UBE3A.iso5,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso5]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso5,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

boxplot(fpkm.stranded.all[UBE3A.iso6,] ~ pheno_stranded$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_stranded)[UBE3A.iso6]), pch=19, xlab="BA4 Motor Cortex", ylab="log2(FPKM+1)", outline=F)
points(fpkm.stranded.all[UBE3A.iso6,] ~ jitter(as.numeric(pheno_stranded$type)), col=as.numeric(pheno_stranded$type))

#Average expression of transcripts
plotMeans("MSTRG.222",bg_stranded, groupvar="type", legend=T)
plotMeans("MSTRG.178",bg_stranded, groupvar="type", legend=T)

dev.off()

#Boxplot comparing expression of the different possible isoform assemblies for antisense UBE3A and UBE3A
pdf(file="ballgown_neurons_output.pdf")

fpkm.neurons.all <- texpr(bg_neurons, meas="FPKM")
fpkm.neurons.all <- log2(fpkm.neurons.all+1)

boxplot(fpkm.neurons.all[iso1.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso1.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso1.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso1.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso1.3,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso1.3,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso1.4,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.4]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso1.4,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso1.5,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso1.5,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[iso2.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso2.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso2.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso2.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso2.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso2.3,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso2.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso2.3,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[iso3.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso3.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso3.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[iso3.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[iso3.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[iso3.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[circRNA,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[circRNA]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[circRNA,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

#UBE3A transcripts
boxplot(fpkm.neurons.all[UBE3A.iso1.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso1.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso1.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso1.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso1.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso1.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[UBE3A.iso2.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso2.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso2.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso2.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso2.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso2.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso2.3,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso2.3]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso2.3,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso2.4,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso2.4]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso2.4,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso2.5,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso2.5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso2.5,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[UBE3A.iso4.1,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso4.1]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso4.1,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))
boxplot(fpkm.neurons.all[UBE3A.iso4.2,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso4.2]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso4.2,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[UBE3A.iso5,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso5]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso5,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

boxplot(fpkm.neurons.all[UBE3A.iso6,] ~ pheno_neurons$type, border=c(1:4), main=paste(ballgown::transcriptNames(bg_neurons)[UBE3A.iso6]), pch=19, xlab="Tissue", ylab="log2(FPKM+1)", outline=F)
points(fpkm.neurons.all[UBE3A.iso6,] ~ jitter(as.numeric(pheno_neurons$type)), col=as.numeric(pheno_neurons$type))

## #Average expression of transcripts
## plotMeans("MSTRG.222",bg_neurons, groupvar="type", legend=T)
## plotMeans("MSTRG.178",bg_neurons, groupvar="type", legend=T)
par(mfrow=c(1,2))
plotTranscripts(ballgown::geneIDs(bg_neurons)[iso1.1], bg_neurons, main=c("MSTRG.178 in Vehicle"),sample=c("vehicle623"))
plotTranscripts(ballgown::geneIDs(bg_neurons)[iso1.1], bg_neurons, main=c("MSTRG.178 in Topotecan"),sample=c("topotecan624"))

plotTranscripts(ballgown::geneIDs(bg_neurons)[UBE3A.iso6], bg_neurons, main=c("MSTRG.222 in Vehicle"),sample=c("vehicle623"))
plotTranscripts(ballgown::geneIDs(bg_neurons)[UBE3A.iso6], bg_neurons, main=c("MSTRG.222 in Topotecan"),sample=c("topotecan624"))

dev.off()
