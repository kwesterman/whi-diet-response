GWAS for putative SFA-genotype interactions influencing LDL-C levels, and construction of a genetic risk score for "responder-ness" to SFA. The analysis is based on the intuition that regression of the product of centered SFA intakes with centered LDL-C levels on allele frequencies can be viewed as a noisy estimation of a latent correlation between SFA and LDL-C. A genome-wide score predicting this correlation would then act as a biomarker of expected SFA response in a given individual.

clean_phenotypes/: Python scripts for assembling phenotypes and covariates.

clean_genotypes/: Bash and plink scripts for processing genotype dosage files into plink2 format for analysis.

run_gwas/: Bash and plink scripts for running the SFA-LDL GWAS.

summarize_gwas.py: Postprocessing of GWAS results (Manhattan/QQ plots, etc.).

meth_interaction.Rmd: Preliminary exploration of a similar concept for DNA methylation (i.e. epigenome-wide association study for the same phenotype).
