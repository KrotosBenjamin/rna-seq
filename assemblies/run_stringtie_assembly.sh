#!/bin/bash
sam_files1='/home/kj/rna-seq/alignments/hg19/unstranded'
sam_files2='/home/kj/rna-seq/alignments/hg19/neurons'
sam_files3='/home/kj/rna-seq/alignments/hg19/stranded'
annotation='/home/kj/rna-seq/assemblies/human/hg19/stranded/strict/annotation.hg19.gtf'
#Assembly with novel annotation file
mkdir -p ballgown/unstranded
mkdir -p ballgown/stranded
mkdir -p ballgown/neurons

stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/cortex432/cortex432.hisat.gtf \
	  $sam_files1/cortex432.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/cortex455/cortex455.hisat.gtf \
	  $sam_files1/cortex455.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/cortex477/cortex477.hisat.gtf \
	  $sam_files1/cortex477.chr15.sorted.bam

stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/heart328/heart328.hisat.gtf \
	  $sam_files1/heart328.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/heart389/heart389.hisat.gtf \
	  $sam_files1/heart389.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/heart435/heart435.hisat.gtf \
	  $sam_files1/heart435.chr15.sorted.bam

stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/liver327/liver327.hisat.gtf \
	  $sam_files1/liver327.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/liver451/liver451.hisat.gtf \
	  $sam_files1/liver451.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/liver463/liver463.hisat.gtf \
	  $sam_files1/liver463.chr15.sorted.bam

stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/lung341/lung341.hisat.gtf \
	  $sam_files1/lung341.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/lung346/lung346.hisat.gtf \
	  $sam_files1/lung346.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/unstranded/lung424/lung424.hisat.gtf \
	  $sam_files1/lung424.chr15.sorted.bam

#Neurons
stringtie -p 8 -G $annotation -e -B -o ballgown/neurons/vehicle623/vehicle623.hisat.gtf \
	  $sam_files2/vehicle623.chr15.sorted.bam

stringtie -p 8 -G $annotation -e -B -o ballgown/neurons/topotecan624/topotecan624.hisat.gtf \
	  $sam_files2/topotecan624.chr15.sorted.bam

#Stranded
stringtie -p 8 -G $annotation -e -B -o ballgown/stranded/ba4.mc830/ba4.mc830.hisat.gtf \
	  $sam_files3/ba4.mc830.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/stranded/ba4.mc831/ba4.mc831.hisat.gtf \
	  $sam_files3/ba4.mc831.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/stranded/ba4.mc832/ba4.mc832.hisat.gtf \
	  $sam_files3/ba4.mc832.chr15.sorted.bam
stringtie -p 8 -G $annotation -e -B -o ballgown/stranded/ba4.mc833/ba4.mc833.hisat.gtf \
	  $sam_files3/ba4.mc833.chr15.sorted.bam
