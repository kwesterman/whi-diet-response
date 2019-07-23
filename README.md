Development of scores predicting response to dietary fat in the Women's Health Initiative. In stage one, a genome-wide interaction study (GWIS) scans for interactions between dietary fat and genotype influencing cardiometabolic risk factors (CRFs). In stage two, these results are aggregated into "fat response scores" and tested in participants in the fat reduction-focused dietary modification trial.

clean_phenotypes/: R and Python scripts for assembling phenotypes and covariates.

clean_genotypes/: Bash and PLINK scripts for processing genotype dosage files into PLINKL format for analysis.

gdi_replication_effort/: Secondary analysis assessing the replication of known GxEs from literature in this dataset.

annotate_res.py: Annotate GWIS results with reference alleles using the associated plinkset.
run_gwis_whi.sh: Script for conducting the genome-wide interaction study using PLINK.
post_gwas.py: Formats GWIS results for downstream analyses and generates Q-Q and Manhattan plots.
pt_model.sh: Generates responder scores from GWIS summary statistics using a series of PLINK calls.
