#!/bin/bash
#De novo assembly
annot='/home/kj/rna-seq/refs/human/hg19/annotation/Genes/genes.gtf'
outputDIR='/home/kj/rna-seq/assemblies/human/hg19/stranded/stict'
mergeLIST='/home/kj/rna-seq/assemblies/human/hg19/stranded/mergelist_gtf.all.txt'

stringtie --merge -p 8 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.0.gtf ./stranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -f 0.05 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.1.gtf ./stranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -F 5 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.2.gtf ./stranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -f 0.05 -F 5 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.3.gtf ./stranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -F 10 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.4.gtf ./stranded/mergelist_gtf.all.txt
stringtie --merge -p 8 -F 5 -T 10 -G $RNA_REFS/human/hg19/annotation/Genes/genes.gtf -o ./stranded/strict/stringtie_merged.5.gtf ./stranded/mergelist_gtf.all.txt
