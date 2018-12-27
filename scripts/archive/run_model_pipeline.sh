#!/bin/bash


POP=$1  # e.g. mesa_white
PLINKSET=$2  # e.g. ../int/mesa


python manhattan.py ../int/sfa/${POP}_res.txt

sbatch --mem 50G -t 1:00:00 pt_model.sh \
	../int/sfa/${POP}_res.txt \
	$PLINKSET/$POP \
	../int/sfa/${POP}_res

sbatch --mem 50G calc_grs.sh \
	../int/sfa/${POP}_res_weights.txt \
	../int/fhs/fhs \
	../int/sfa/fhs_scores_${POP}_model

