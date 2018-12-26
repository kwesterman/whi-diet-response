#!/bin/bash

module load plink

WEIGHTSFILE=$1
PLINKSET=$2
OUTFILE=$3

plink --bfile $PLINKSET \
	--score $WEIGHTSFILE \
	--out $OUTFILE
