#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunDataPreparations"

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load java-1.8.0_40
module load python/2.7.15/native
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=[Path to tool folder with nextflow]

NXF_VER=20.10.0 ${nextflow_path}/nextflow run PerCohortDataPreparations.nf \
--genopath '[Folder with genotype files in .h5 format]' \
--expressionpath '[Prepared and normalised gene expression file]' \
--covariatepath '[Covariate file]' \
--gte '[File with sample ID matches between genotype and expression data]' \
--outputpath '[Folder where to write encoded files]' \
--studyname '[CohortName_GeneExpressionPlatform]' \
--numcovariates [Nr of covariates to include] \
--profile [slurm,singularity] \
-with-report PerCohortDataPreparationsReport.html \
-resume
