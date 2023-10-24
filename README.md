# Per cohort data preparations

This is Nextflow pipeline per-cohort data preparations for running encoded  genome-wide eQTL meta-analyses.

Specifically, this pipeline performs the following steps:

- Uses 1000G 30x reference file to create variant mapper files (to make the SNPs from different studies jointly analyzable in central site).
- Encodes genotype and gene expression data, and deletes the random matrix used for encoding. This means that no information for individual study participant is obtainable from encoded data, even for the original cohort analyst.
- Calculates partial derivatives, needed for running the eQTL mapping.
- Permutes the sample links on the unencoded data, encodes and calculates partial derivatives for the permuted data. This is needed for obtaining in-sample LD estimates for downstream analyses in the central site (useful for e.g. multiple testing corrections and fine-mapping).
- Associates expression PCs with genotypes, writes out suggestive associations (P<1×10^-5^). This enables us to make informed decision which covariates to include into encoded HASE model in the central site and control to for the collider effects.
- Replaces original sample IDs in the encoded data with "CohortName_index".
- Collects several summary reports, QC reports and diagnostic plots from the output of [data QC pipeline](#1-data-qc).
- Organises all the partial derivates, encoded matrices and QC reports into the structured folder structure, ready for sharing.
- Calculates `md5sum` for all the shared files, so the integrity of the upload can be checked.

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java 8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline
You can either clone it by using git (if available in HPC):
`git clone https://github.com/eQTLGen/PerCohortDataPreparatons.git`

Or just download this from the gitlab/github download link and unzip.

### Input files

- Imputed genotype data in `.hdf5` format. You can convert your imputed `.vcf.gz` files into `.hdf5` format by `ConvertVcf2Hdf5` pipeline.
- Standardized folder which contains gene expression matrix, covariate matrix and data QC summary files. This is the output of `dataqc` pipeline, see the details from the documentation of [Data QC pipeline](TBA).

### Running the data preparation command

Go to folder `PerCohortPreparations` and modify the script template `submit_per_cohort_preparations_template.sh` with your input paths. Below is an example template for Slurm scheduler. Some of the paths are pre-filled, assuming that you follow [eQTLGen phase II cookbook](https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook) and its recommended folder structure, however you can also use custom paths.

```bash
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
```

You can save the modified script version to informative name, e.g. `submit_per_cohort_preparations_[**CohortName_PlatformName**].sh`.

You can select HPC scheduler type by adjusting the profile as following:

- Slurm: `-profile slurm,singularity`
- PBS/TORQUE: `-profile pbs,singularity`
- SGE: `-profile sge,singularity`

Then submit the job `sbatch submit_per_cohort_preparations_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.

### Monitoring and debugging

- Monitoring:
  - Monitor the `slurm-***.out` log file and check if all the steps finish without error. Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.
  - Use `squeue -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file  is automatically written to your output folder `pipeline_info` subfolder, for potential errors or warnings. E.g. `output/pipeline_info/DataQcReport.html`.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from `work` folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted.
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.

### Output

- In the `output` folder there is subfolder called `[YourCohortName]_IntermediateFilesEncoded_to_upload`. This folder contains all the non-personal files, logs, reports and should be shared with the central site.  This includes: encoded phenotype data (encoded gene expression matrix), encoded genotype data (encoded genotype files), partial derivatives folder, folder containing SNP QC information file, encoded files for permuted data, partial derivatives for permuted data, several diagnostic plots and summary statistic files for gene expression summaries.
- In the `output` folder there is also file `[YourCohortName]_IntermediateFilesEncoded_to_upload.md5`. This can be used for checking the integrity of the files in folder.

## Acknowledgements

This pipeline utilizes HASE (https://github.com/roshchupkin/hase) and some of its helper scripts originally developed by:

Gennady V. Roscupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands) 

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands). 

### Changes from the original HASE repo

Robert Warmerdam (Department of Genetics, University Medical Center Groningen, University of Groningen, Groningen, Netherlands) modified the original HASE and fixed some bugs.

Urmo Võsa (Institute of Genomics, University of Tartu, Tartu, Estonia) incorporated it into Nextflow pipeline and customized/supplanted some parts of the code.

**Changes:**

- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the `--intercept` option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Setting /hdgwas/data.py were changed `self.chunk_size=10000 --> self.chunk_size=20000`.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep the tool lightweight.
- In HASE "data" folder, .txt file with link to original 1000G reference was removed, to avoid confusion in eQTLGen phase II analysis.

### Citation

Original method paper for HASE framework:

[Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)](https://www.nature.com/articles/srep36076)

### Contacts

For this Nextflow pipeline: urmo.vosa at gmail.com

For the method of HASE, find contacts from [original HASE repo](https://github.com/roshchupkin/hase)

