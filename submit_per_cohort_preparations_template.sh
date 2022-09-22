#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunDataPreparations"

# These are needed modules in UT HPC to get Singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

nextflow_path=../../tools

genotypes_hdf5=../../3_ConvertVcf2Hdf5/output # Folder with genotype files in .hdf5 format
qc_data_folder=../../1_DataQC/output # Folder containing QCd data, inc. expression and covariates
output_path=../output

NXF_VER=21.10.6 ${nextflow_path}/nextflow run PerCohortDataPreparations.nf \
--hdf5 ${genotypes_hdf5} \
--qcdata ${qc_data_folder} \
--outdir ${output_path} \
-profile slurm,singularity \
-resume
