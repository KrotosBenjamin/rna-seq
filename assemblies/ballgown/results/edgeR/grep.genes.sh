#!/bin/bash

# grep -Fwf genes.txt exon.DE_unstranded_cortexVSheart.txt | \
#     sed -r "s/\s+/\t/g" > exon_DE_unstranded_cortexVSheart.tsv
grep -Fwf genes.txt ./unstranded/exon.DE_unstranded_cortexVSliver.txt | \
    sed -r "s/\s+/\t/g" > ./unstranded/exon_DE_unstranded_cortexVSliver.tsv
grep -Fwf genes.txt ./unstranded/exon.DE_unstranded_cortexVSlung.txt | \
    sed -r "s/\s+/\t/g" > ./unstranded/exon_DE_unstranded_cortexVSlung.tsv
