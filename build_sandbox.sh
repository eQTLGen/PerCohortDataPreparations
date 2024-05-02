#!/bin/bash

#SBATCH --time=00:59:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="Prepare singularity sandbox container"

# These are needed modules in UT HPC to get Singularity and Nextflow running.
# Replace with appropriate ones for your HPC.
module load singularity/3.11.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR="../../singularitycache"
export NXF_HOME="../../nextflowcache"
export TMPDIR="../../singularitytmp"
export SINGULARITY_TMPDIR="../../singularitytmp"

singularity build --sandbox singularity_img/quay.io-eqtlgen-eqtlgen_phase2-v1.8-sbox singularity_img/quay.io-eqtlgen-eqtlgen_phase2-v1.8.img