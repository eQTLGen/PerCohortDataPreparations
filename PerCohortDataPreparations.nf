#!/usr/bin/env nextflow

// Things this pipeline does:
// 1. Applies inverse normal transformation on each gene in the input file.
// 2. Replaces the sample IDs in the input gene expression file with those as provided in genotype-to-expression file.
// 3. Replaces the probe names with gene names from empirical probe matching approach.
// 4. Chunks gene expression file into chunks.
// 5. Constructs mapper files for genotype data (currently assumes the file named explicitly 1000Gp1v3.ref.gz in the data folder of tool path).
// 6. Encodes the data.
// 7. Calculates partial derivatives.
// 8. Replaces sample IDs in the encoded data with mock IDs in the form of Cohort__SampleNumber.
// 9. Organise files needed for sharing.


def helpmessage() {

log.info"""

PerCohortDataPreparations
==============================================
Pipeline for running data preparation and encoding steps in individual datasets, in order to enable centralized meta analysis over all datasets.

This pipeline mostly recapitulates the steps as outlined in HASE wiki page (https://github.com/roshchupkin/hase/wiki/HD-Meta-analysis-GWAS), with some bugfixes.

Usage:

nextflow run PerCohortDataPreparations.nf \
--genopath '/Genotype folder/' \
--expressionpath '[Expression file]' \
--covariatepath '[Covariate file]' \
--probematches '[File with array probe matches]' \
--gte '[Genotype to expression file]' \
--studyname '[NameOfYourStudy]' \
--outputpath '[/Output folder/]'\
--numcovariates 20

Mandatory arguments:
--genopath        Path to input genotype folder. It has to be in hdf5 format and contain at least 3 subfolders (genotypes, probes, individuals).
--expressionpath  Path to the gene expression file. May be gzipped.
--covariatepath   Path to covariate file. May be gzipped.
--probematches    Path to the file connecting array probes with respective ENSEMBL gene IDs.
--gte             Path to the file connecting genotype data sample IDs with gene expression matrix sample IDs.
--studyname       Name of the study.
--outputpath      Path to output folder where encoded, prepared and organised data is written.
--numcovariates   Number of covariates to correct the analysis for. It defaults to 20 first columns in the covariate file.


""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.studyname = ''
params.genopath = ''
params.numcovariates = ''
params.expressionpath = ''
params.covariatepath = ''
params.probematches = ''
params.gte = ''
params.outputpath = ''

log.info """================================================================
HASE per-cohort preparation pipeline v${workflow.manifest.version}"
================================================================"""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Input genotype directory']                 = params.genopath
summary['Input gene expression file']               = params.expressionpath
summary['Input covariate file']                     = params.covariatepath
summary['Input probe mapping file']                 = params.probematches
summary['Input genotype-to-expression file']        = params.gte
summary['Number of covariates']                     = params.numcovariates
summary['Output directory']                         = params.outputpath
summary['Study names']                              = params.studyname

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

Genotypes = Channel.fromPath(params.genopath)
Genotypes.into{genotypes_to_mapper; genotypes_to_encoding; genotypes_to_pd; genotypes_to_genpc; genotypes_to_organise}

expression = Channel.fromPath(params.expressionpath)
covariates = Channel.fromPath(params.covariatepath)

process PrepareExpressionData {

    input:
      path expression from expression
      path emp from params.probematches
      path gte from params.gte

    output:
      path './exp_data/' into expression_to_encoding, expression_to_pd

    """
    python $baseDir/bin/helperscripts/convert_measurement_data.py \
    -i ${expression} \
    -gte ${gte} \
    -o ./exp_data/ \
    -n 20000
    """
}

process PrepareCovariateData {

    input:
      path covariates from covariates
      path emp from params.probematches
      path gte from params.gte
      val num_of_covariates from params.numcovariates

    output:
      path './covariate_data/' into covariates_to_pd

    """
    python $baseDir/bin/helperscripts/convert_measurement_data.py \
    -i ${covariates} \
    -gte ${gte}\
    -o ./covariate_data/ \
    -n 1 \
    -nc ${num_of_covariates}
    """

}

process CreateMapperFiles {

    input:
      path genopath from genotypes_to_mapper
      val studyname from params.studyname

    output:
      path './mapper/' into mapper_to_encode, mapper_to_pd, mapper_to_organize

    """
    python $baseDir/bin/hase/tools/mapper.py \
    -g ${genopath} \
    -study_name ${studyname} \
    -o ./mapper/ \
    -chunk 35000000 \
    -ref_name 1000Gp1v3_ref \
    -ref_chunk 1000000 \
    -probe_chunk 1000000
    """
}

process EncodeData {

    input:
      path mapper from mapper_to_encode
      path genopath from genotypes_to_encoding
      path expression from expression_to_encoding
      val studyname from params.studyname

    output:
      file './encoded/' into encoded

    """
    echo ${genopath}
    echo ${mapper}
    echo ${expression}
    echo ${params.studyname}

    python $baseDir/bin/hase/hase.py \
    -g ${genopath} \
    -study_name ${studyname} \
    -o ./encoded/ \
    -mapper ${mapper}/ \
    -ph ${expression} \
    -mode encoding

    # Remove random matrices to make back-encoding impossible
    rm /encoded/F*
    """
}

process PartialDerivatives {

    input:
      path mapper from mapper_to_pd
      path genopath from genotypes_to_pd
      path expression from expression_to_pd
      path covariates from covariates_to_pd
      val studyname from params.studyname

    output:
      file './pd/' into pd

    """
    python $baseDir/bin/hase/hase.py \
    -g ${genopath}/ \
    -study_name ${studyname} \
    -ph ${expression}/ \
    -cov ${covariates}/ \
    -mapper ${mapper}/ \
    -o ./pd/ \
    -mode single-meta
    """
}

process PrepareGenRegPcs {
    // TODO: continue here

    input:
      path covariates from covariates_to_pd

    output:
      path cov_folder
      path pheno_folder

    """
    # Split covariate file to two pieces: covariates (4 first MDS) and 100 first PCs
    awk -F'\t' '{ print $1, $2, $3, $4 }' > covariate_MDS.txt
    awk -F'\t' '{ print $1, $5..$105 }' > pheno_expPC.txt

    mkdir cov_folder
    mkdir pheno_folder

    mv covariate_MDS.txt cov_folder/.
    mv pheno_expPC.txt pheno_folder/.
    """
}

process RunGenRegPcs {

    input:
      path genopath from genotypes_to_genpc
      path covariates from cov_folder
      path phenotypes from pheno_folder
      val studyname from params.studyname
      each Chunk from 1..100

    output:
      *_GenRegPcs.txt into genetic_pcs

    """
    python $baseDir/bin/hase/hase.py \
    -g ${genopath} \
    -study_name ${studyname} \
    -o output \
    -ph pheno_folder \
    -cov cov_folder \
    -th 4.42 \  # for now hardcoded threshold, corresponding to ~2*10-5 for N=100 and 9.974149e-06 for N=10000
    -mode regression \ 
    -maf 0.0 \
    -node 100 ${Chunk} \ # where N is number of current job, and 100 is total number of jobs. 
    -cluster "y"

    python $baseDir/bin/hase/tools/analyzer.py \ 
    -r output \
    -o ${Chunk}_GenRegPcs.txt
    """
}

process OrganizeEncodedData {

    input:
      path pd from pd
      path encoded from encoded
      path mapper from mapper_to_organize
      path genopath from genotypes_to_organise
      val studyname from params.studyname
      path genetic_pcs from genetic_pcs.collectFile(name: 'GenRegPcs.txt', keepHeader: true, sort: true)

    output:
      path '*' into OrganizedFiles

    """
    echo ${pd}
    echo ${encoded}
    echo ${mapper}

    mkdir -p IntermediateFilesEncoded_to_upload/EncodedGenotypeData/genotype
    mkdir IntermediateFilesEncoded_to_upload/EncodedGenotypeData/individuals
    mkdir IntermediateFilesEncoded_to_upload/EncodedPhenotypeData
    mkdir IntermediateFilesEncoded_to_upload/pd_shared

    cp -r ${genopath}/probes IntermediateFilesEncoded_to_upload/EncodedGenotypeData/
    cp -r ${genopath}/SNPQC IntermediateFilesEncoded_to_upload/

    cp -r ./${encoded}/encode_individuals/* IntermediateFilesEncoded_to_upload/EncodedGenotypeData/individuals/
    cp -r ./${encoded}/encode_genotype/* IntermediateFilesEncoded_to_upload/EncodedGenotypeData/genotype/
    cp -r ./${encoded}/encode_phenotype/* IntermediateFilesEncoded_to_upload/EncodedPhenotypeData/
    cp -r ./${mapper} IntermediateFilesEncoded_to_upload/

    cp ./${pd}/*.npy IntermediateFilesEncoded_to_upload/pd_shared/

    mv GenRegPcs.txt IntermediateFilesEncoded_to_upload/{studyname}_GenRegPcs.txt
    """
}

process ReplaceSampleNames {

    publishDir "${params.outputpath}", mode: 'copy', overwrite: true

    input:
      path OrganizedFiles from OrganizedFiles
      val studyname from params.studyname

    output:
      path IntermediateFilesEncoded_to_upload into IntermediateFilesEncodedSampleIdsReplaced_to_upload

    """
    python $baseDir/bin/helperscripts/replace_sample_names.py -IntFileEnc ${OrganizedFiles}

    # Calculate md5sum for all the files to share and add write this out
    find -type f -exec md5sum '{}' \; ; > ${studyname}_OrganizedFiles.md5
  
    """
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}