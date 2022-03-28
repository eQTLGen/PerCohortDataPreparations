#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunDataPreparations"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=../../tools

genotypes_hdf5=[Folder with genotype files in .hdf5 format]
qc_data_folder=[Folder containing QCd data, inc. expression and covariates]
output_path=../output

NXF_VER=21.10.6 ${nextflow_path}/nextflow run PerCohortDataPreparations.nf \
--hdf5 ${genotypes_hdf5} \
--qcdata ${qc_data_folder} \
--outdir ${output_path} \
-profile slurm,singularity \
-resume
