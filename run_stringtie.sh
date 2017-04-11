#!/bin/bash
#De novo assembly
sam_files='/home/kj/rna-seq/alignments/hg19'
mkdir unstranded

stringtie -p 8 -l cortex432 -o unstranded/cortex432/transcripts.chr15.gtf \
	  $sam_files/cortex432.chr15.sorted.bam
stringtie -p 8 -l cortex455 -o unstranded/cortex455/transcripts.chr15.gtf \
	  $sam_files/cortex455.chr15.sorted.bam
stringtie -p 8 -l cortex477 -o unstranded/cortex477/transcripts.chr15.gtf \
	  $sam_files/cortex477.chr15.sorted.bam

stringtie -p 8 -l heart328 -o unstranded/heart328/transcripts.chr15.gtf \
	  $sam_files/heart328.chr15.sorted.bam
stringtie -p 8 -l heart389 -o unstranded/heart389/transcripts.chr15.gtf \
	  $sam_files/heart389.chr15.sorted.bam
stringtie -p 8 -l heart435 -o unstranded/heart435/transcripts.chr15.gtf \
	  $sam_files/heart435.chr15.sorted.bam

stringtie -p 8 -l liver327 -o unstranded/liver327/transcripts.chr15.gtf \
	  $sam_files/liver327.chr15.sorted.bam
stringtie -p 8 -l liver451 -o unstranded/liver451/transcripts.chr15.gtf \
	  $sam_files/liver451.chr15.sorted.bam
stringtie -p 8 -l liver463 -o unstranded/liver463/transcripts.chr15.gtf \
	  $sam_files/liver463.chr15.sorted.bam

stringtie -p 8 -l lung341 -o unstranded/lung341/transcripts.chr15.gtf \
	  $sam_files/lung341.chr15.sorted.bam
stringtie -p 8 -l lung346 -o unstranded/lung346/transcripts.chr15.gtf \
	  $sam_files/lung346.chr15.sorted.bam
stringtie -p 8 -l lung424 -o unstranded/lung424/transcripts.chr15.gtf \
	  $sam_files/lung424.chr15.sorted.bam

#Merge files
ls -1 unstranded/*/transcripts.chr15.gtf > unstranded/mergelist_gtf.all.txt
ls -1 unstranded/*rt*/transcripts.chr15.gtf > unstranded/mergelist_gtf.hVSc.txt
ls -1 unstranded/cortex*/transcripts.chr15.gtf > unstranded/mergelist_gtf.cortex.txt

stringtie --merge -p 8 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf \
	  -o unstranded/stringtie_merged.all.gtf unstranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf \
	  -o unstranded/stringtie_merged.hVSc.gtf unstranded/mergelist_gtf.hVSc.txt
stringtie --merge -p 8 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf \
	  -o unstranded/stringtie_merged.cortex.gtf unstranded/mergelist_gtf.cortex.txt

# #Second assembly
# mkdir unstranded/ballgown

# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/cortex432/cortex432.hisat.gtf \
# 	  $sam_files/cortex432.chr15.sorted.bam
# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/cortex455/cortex455.hisat.gtf \
# 	  $sam_files/cortex455.chr15.sorted.bam
# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/cortex477/cortex477.hisat.gtf \
# 	  $sam_files/cortex477.chr15.sorted.bam
# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/heart328/heart328.hisat.gtf \
# 	  $sam_files/heart328.chr15.sorted.bam
# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/heart389/heart389.hisat.gtf \
# 	  $sam_files/heart389.chr15.sorted.bam
# stringtie -p 8 -G unstranded/stringtie_merged.gtf -e -B \
# 	  -o unstranded/ballgown/heart435/heart435.hisat.gtf \
# 	  $sam_files/heart435.chr15.sorted.bam
