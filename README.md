# Per cohort data preparations

This is Nextflow pipeline per-cohort data preparations for running encoded pooled genome-wide eQTL analyses.

## Usage instructions

### Requirements for the system

- Have access to HPC with multiple cores.

- Have Bash >=3.2 installed.

- Have Java 8 installed.

- Have SLURM scheduler managing the jobs in the HPC.

- Have Singularity (preferred) and/or conda installed.

### Requirements for the input data

- Imputed genotype data: This is the output of the ConvertVcf2Hdf5 pipeline and should be in hdf5 format.

- Expression file:

	- Tab-delimited expression matrix.

	- Each row contains gene expression value and each column contains sample ID.

	- First column has header `-` and contains either ENSEMBL gene IDs (RNA-seq), Illumina ArrayAddress IDs for probes or probe names for corresponding Affymetrix array.

	- File can be gzipped.

	- Expression QC has been done.

	- Data has been appropriately normalised and transformed.

- Covariate file:

	- Tab delimited matrix.

	- Each row contains sample ID and each column contains covariate.

	- First column has header `-` and contains sample ID.

	- File can be gzipped.

	- In the pipeline you can specify how many covariates to include into data preparation.

- Probe matches file: tab-delimited file with matches between array probes and corresponding ENSEMBL gene IDs. Column names are `Probe` and `Ensembl`. If you already have ENSEMBL IDs in your data, then you can just use file with two columns of ENSEMBL IDs. For eQTLGen phase 2 analyses we will provide you the needed empirical probe matching files for each platform.

- Genotype-to-expression file: tab-delimited file with sample ID matches between genotype data (column 1) and expression data (column 2). If those IDs are the same in both data layers, you can just use file with same sample IDs in both columns.

### Setup the pipeline

1. Make a folder named `per_cohort_data_preparations` and get in there.

2. Make subfolder `tools`, go into this folder and have Nextflow executable installed: https://www.nextflow.io/docs/latest/getstarted.html#installation. After that, step back up into `per_cohort_data_preparations`.

3. Get the genotype conversion pipeline from here: https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5

You can either clone it by using git (if available in HPC):

`git clone https://gitlab.com/eqtlgen-group/ConvertVcf2Hdf5`

Or just download this from the gitlab download link and unzip into `per_cohort_data_preparations` folder.

4. Put genotype reference file into `PerCohortPreparations/bin/hase/data/` folder. For eQTLGen phase 2 analyses, you can get this file from **here**.

5. Make folder where you want to save your encoded files for sharing. This can be in the `per_cohort_data_preparations` folder or elsewhere.

### Running the conversion command

Go to folder `PerCohortPreparations` and modify the script template `submit_per_cohort_preparations_template.sh` with your input paths and .

```
#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --job-name="RunHASEwNextflow"

# Load needed system tools (Java 8 is required, one of singularity or anaconda - python 2.7 is needed,
# depending on the method for dependancy management). The exact names of tool modules might depend on HPC.
module load java-1.8.0_40
module load python/2.7.15/native
module load singularity/3.5.3
module load squashfs/4.4

nextflow_path=../tools/

NXF_VER=20.10.0 ${nextflow_path}/nextflow run PerCohortDataPreparations.nf \
--genopath '[Folder with genotype files in .h5 format]' \
--expressionpath '[Prepared and normalised gene expression file]' \
--covariatepath '[Covariate file]' \
--probematches '[File with matches between array probe IDs and gene names]' \
--gte '[File with sample ID matches between genotype and expression data]' \
--outputpath '[Folder where to write encoded files]' \
--studyname '[CohortName_GeneExpressionPlatform]' \
--NrOfCovariates [Nr of covariates to include] \
--profile '[cluster_slurm,singularity_profile/conda_profile]' \
-with-report PerCohortPreparationsReport.html \
-resume

```

You can save the modified script version to informative name, e.g. `submit_per_cohort_preparations_[**CohortName_PlatformName**].sh`.

Then submit the job `sbatch submit_per_cohort_preparations_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.

### Monitoring and cleaning up

- Monitor the `slurm-***.out` log file and check if all steps finish without error.

	Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.

- If pipeline crashes and you figure out how to fix this, you can just resubmit the same script after fixes. Nextflow does not rerun unaffected steps and continues only from the steps affected by the issue. 

- When the work has finished, check `PerCohortPreparationsReport.html` for potential errors or warnings.

- For some analyses, `work` directory can become quite large. So, when the pipeline is successfully finished, you can delete it in order to save disk space of your HPC. But this means that the pipeline will restart from scratch, if you ever need to rerun it.

### Results of the pipeline

After successful completion of the pipeline, there should be subfolder named `IntermediateFilesEncoded_to_upload` under your output folder. This includes: encoded phenotype data (encoded gene expression matrix), encoded genotype data (encoded genotype files), partial derivatives folder, and folder containing SNP QC information file.

### Benchmark

Here is the estimate, how much time the conversion is expected to take.

- Initial number of samples in .vcf genotype data: \~10,000

- Number of samples to filter into .h5 genotype format: \~1,000

- Infrastructure: University HPC with \~150 compute nodes

- Dependency management: Singularity 

- Time to run the pipeline (without wall times): \~4.5h

- CPU hours: \~80h

- Final size of `work` subdirectory: \~500GB\

## Acknowledgements

This pipeline utilizes HASE (https://github.com/roshchupkin/hase) and some of its helper scripts originally developed by:

Gennady V. Roscupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands) 

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands). 

### Changes from the original HASE repo

Robert Warmerdam (Department of Genetics, University Medical Center Groningen, University of Groningen, Groningen, Netherlands) modified the original HASE and fixed some bugs.

Urmo Võsa (Institute of Genomics, University of Tartu, Tartu, Estonia) incorporated it into Nextflow pipeline and customized some parts of the code.

**Changes:**

- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the `--intercept` option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Lines 355-357 of hase.py were commented out because this caused pipeline to crash when >1 datasets were added.
- Line 355 of /hdgwas/data.py were changed `self.chunk_size=10000 --> self.chunk_size=20000`.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep the tool lightweight.

### Citation

Original method paper for HASE framework:

[Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)](https://www.nature.com/articles/srep36076)

### Contacts

For this Nextflow pipeline: urmo.vosa@gmail.com

For the method of HASE, find contacts from here: https://github.com/roshchupkin/hase
