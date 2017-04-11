#!/bin/bash
#De novo assembly
sam_files='/home/kj/rna-seq/alignments/hg19/stranded'

stringtie -p 8 -l ba4.mc830 -o stranded/ba4.mc830/transcripts.chr15.gtf \
	  $sam_files/ba4.mc830.chr15.sorted.bam
stringtie -p 8 -l ba4.mc831 -o stranded/ba4.mc831/transcripts.chr15.gtf \
	  $sam_files/ba4.mc831.chr15.sorted.bam
stringtie -p 8 -l ba4.mc832 -o stranded/ba4.mc832/transcripts.chr15.gtf \
	  $sam_files/ba4.mc832.chr15.sorted.bam
stringtie -p 8 -l ba4.mc833 -o stranded/ba4.mc833/transcripts.chr15.gtf \
	  $sam_files/ba4.mc833.chr15.sorted.bam

#Merge files
ls -1 stranded/*/transcripts.chr15.gtf > stranded/mergelist_gtf.all.txt

# stringtie --merge -p 8 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf \
# 	  -o stranded/stringtie_merged.gtf stranded/mergelist_gtf.txt

# #Second assembly
# mkdir stranded/ballgown

# stringtie -p 8 -G stranded/stringtie_merged.gtf -e -B \
# 	  -o stranded/ballgown/ba4.mc830/ba4.mc830.hisat.gtf \
# 	  $sam_files/ba4.mc830.chr15.sorted.bam
# stringtie -p 8 -G stranded/stringtie_merged.gtf -e -B \
# 	  -o stranded/ballgown/ba4.mc831/ba4.mc831.hisat.gtf \
# 	  $sam_files/ba4.mc831.chr15.sorted.bam
# stringtie -p 8 -G stranded/stringtie_merged.gtf -e -B \
# 	  -o stranded/ballgown/ba4.mc832/ba4.mc832.hisat.gtf \
# 	  $sam_files/ba4.mc832.chr15.sorted.bam
# stringtie -p 8 -G stranded/stringtie_merged.gtf -e -B \
# 	  -o stranded/ballgown/ba4.mc833/ba4.mc833.hisat.gtf \
# 	  $sam_files/ba4.mc833.chr15.sorted.bam
