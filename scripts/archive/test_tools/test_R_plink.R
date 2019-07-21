library(foreach)
library(itertools)
library(broom)

Rplink <- function(PHENO, GENO, CLUSTER, COVAR) {

	f1 <- function(v) {
		fit <- tidy(lm(rnorm(PHENO ~ v))
		m <- c(fit$estimate[2], fit$p.value[2])
		m <- head(PHENO)
		c(length(m), m)
	}

	foreach(g=iter(GENO, by="col"), .combine=cbind) %do% f1(g)
}

# Run with: plink --bfile mytest --pheno /cluster/home/kweste01/kw/diet_response/data/processed/gen3/whi_white_gwas_phenos.txt --mpheno 14 --R test_R_plink.R --out mytest
# Seems to require phenotype to be specifically provided/not just a set of zeros (?)
# Must return a numeric vector
# Still haven't really seen this produce any actual numeric results in the resulting *.auto.R file?
