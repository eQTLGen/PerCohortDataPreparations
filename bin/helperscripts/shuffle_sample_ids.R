args <- commandArgs(trailingOnly = TRUE)

library(data.table)

# do not set seed
# set.seed(123)

exp <- fread(args[1], header = TRUE)
cov <- fread(args[2], header = TRUE)

# Take the intersection of samples over two datasets
exp <- exp[exp$ID %in% cov$ID, ]
cov <- cov[cov$ID %in% exp$ID, ]
# Order the datasets based on sample ID, so that link between phenotype and covariate sample IDs remains 
# after shuffling
exp <- exp[order(exp$ID), ]
cov <- cov[order(cov$ID), ]

if(!all(exp$ID == cov$ID)){stop("Error: phenotype and covariate sample IDs are not aligned!")}

exp_suff <- exp

# For covariates: keep the link between expression IDs and expression PC IDs and shuffle the link with genotypes.
# For 10 first MDSs keep the link with genotypes unchanged and shuffle the link with expression IDs.
cov_suff_expPC <- cov[, c(1, 11:111), with = FALSE]
cov_suff_genPC <- cov[, c(1, 2:11), with = FALSE]

exp_suff$ID <- sample(exp_suff$ID, length(exp_suff$ID))
cov_suff_expPC$ID <- exp_suff$ID
cov_suff <- merge(cov_suff_genPC, cov_suff_expPC, by = "ID")

fwrite(exp_suff, args[3], sep = "\t", row.names = FALSE)
fwrite(cov_suff, args[4], sep = "\t", row.names = FALSE)
