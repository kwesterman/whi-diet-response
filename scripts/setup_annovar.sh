#!/bin/bash


PREFIX=$1  
INPUTFILE=$PREFIX.weights  # *.weights file from pipeline
OUTFILE=$PREFIX.avinput

tail -n +2 $INPUTFILE | awk -v OFS='\t' '{print $13,$14,$14,$4,$5, $1,$9}' > $OUTFILE
