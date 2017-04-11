#!/bin/bash

grep -Fwf transcripts.txt DE_unstranded_cortexVSliver.csv | \
    sed -r "s/\s+/\t/g" > DE_unstranded_cortexVSliver.tsv
grep -Fwf transcripts.txt DE_unstranded_cortexVSlung.csv | \
    sed -r "s/\s+/\t/g" > DE_unstranded_cortexVSlung.tsv
grep -Fwf transcripts.txt DE_unstranded_heartVScortex.csv | \
    sed -r "s/\s+/\t/g" > DE_unstranded_heartVScortex.tsv
