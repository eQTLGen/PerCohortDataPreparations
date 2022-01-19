#!/usr/bin/env nextflow

// Things this pipeline does:
// 1. Replaces the sample IDs in the input gene expression file with those as provided in genotype-to-expression file.
// 2. Chunks gene expression file into chunks.
// 3. Constructs mapper files for genotype data (currently assumes the file named explicitly 1000Gp1v3.ref.gz in the data folder of tool path).
// 4. Encodes the data.
// 5. Calculates partial derivatives.
// 6. Replaces sample IDs in the encoded data with mock IDs in the form of Cohort__SampleNumber.
// 7. Organise files needed for sharing.


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
--gte '[Genotype to expression file]' \
--studyname '[NameOfYourStudy]' \
--outputpath '[/Output folder/]'\
--numcovariates 20

Mandatory arguments:
--genopath        Path to input genotype folder. It has to be in hdf5 format and contain at least 3 subfolders (genotypes, probes, individuals).
--expressionpath  Path to the gene expression file. May be gzipped.
--covariatepath   Path to covariate file. May be gzipped.
--gte             Path to the file connecting genotype data sample IDs with gene expression matrix sample IDs.
--studyname       Name of the study. Needs to match with the study name as in the file names in genotype .hdf5 folder.
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
summary['Input genotype-to-expression file']        = params.gte
summary['Number of covariates']                     = params.numcovariates
summary['Output directory']                         = params.outputpath
summary['Study names']                              = params.studyname

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

Genotypes = Channel.fromPath(params.genopath)
Genotypes.into{genotypes_to_mapper; genotypes_to_encoding; genotypes_to_pd; genotypes_to_perm_encoding; genotypes_to_perm_pd; genotypes_to_genpc; genotypes_to_organise}

expression = Channel.fromPath(params.expressionpath)
covariates = Channel.fromPath(params.covariatepath)

process PrepareExpressionData {

    input:
      path expression from expression
      path gte from params.gte

    output:
      path './exp_data/' into expression_to_encoding, expression_to_pd, expression_to_permutation

    """
    mkdir exp_data

    python $baseDir/bin/helperscripts/convert_measurement_data.py \
    -t exp \
    -i ${expression} \
    -gte ${gte} \
    -o ./exp_data/ \
    -n 20000
    """
}

process PrepareCovariateData {

    input:
      path covariates from covariates
      path gte from params.gte
      val num_of_covariates from params.numcovariates

    output:
      path './covariate_data/' into covariates_to_pd, covariates_to_permutation, covariates_to_genpc

    """
    mkdir covariate_data

    python $baseDir/bin/helperscripts/convert_measurement_data.py \
    -t cov \
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
      path './mapper/' into mapper_to_encode, mapper_to_pd, mapper_to_perm_encode, mapper_to_perm_pd, mapper_to_organize

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
    rm ./encoded/F*
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

process PermuteData {

  input:
    path expression from expression_to_permutation
    path covariates from covariates_to_permutation

  output:
    path 'shuffled_expression_folder' into perm_exp_to_encoding, perm_exp_to_pd
    path 'shuffled_covariates_folder' into perm_cov_to_encoding, perm_cov_to_pd

  """
  mkdir shuffled_expression_folder
  mkdir shuffled_covariates_folder

  Rscript --vanilla $baseDir/bin/helperscripts/shuffle_sample_ids.R \
  ${expression}/batch0.txt \
  ${covariates}/covariates_HASE_format.txt \
  ./shuffled_expression_folder/shuffled_expression.txt \
  ./shuffled_covariates_folder/shuffled_covariates.txt
  """
}

process EncodeDataPermuted {

    input:
      path mapper from mapper_to_perm_encode
      path genopath from genotypes_to_perm_encoding
      path expression from perm_exp_to_encoding
      val studyname from params.studyname

    output:
      file './encoded_permuted/' into encoded_permuted

    """
    python $baseDir/bin/hase/hase.py \
    -g ${genopath} \
    -study_name ${studyname} \
    -o ./encoded/ \
    -mapper ${mapper}/ \
    -ph ${expression} \
    -mode encoding

    # Remove random matrices to make back-encoding impossible
    rm ./encoded/F*

    mv encoded encoded_permuted
    """
}


process PartialDerivativesPermuted {

    input:
      path mapper from mapper_to_perm_pd
      path genopath from genotypes_to_perm_pd
      path expression from perm_exp_to_pd
      path covariates from perm_cov_to_pd
      val studyname from params.studyname

    output:
      file './pd_permuted/' into pd_permuted

    """
    python $baseDir/bin/hase/hase.py \
    -g ${genopath}/ \
    -study_name ${studyname} \
    -ph ${expression}/ \
    -cov ${covariates}/ \
    -mapper ${mapper}/ \
    -o ./pd/ \
    -mode single-meta

    mv pd pd_permuted
    """
}

process PrepareGenRegPcs {

    input:
      path covariates from covariates_to_genpc

    output:
      path cov_folder
      path pheno_folder

    """
    # Split covariate file to two pieces: covariates (10 first MDS) and 100 first PCs
    awk -F'\t' '{ print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11}' ${covariates}/covariates_HASE_format.txt > covariate_MDS.txt
    awk 'BEGIN{FS=OFS="\t"}{printf \$1"\t"}{for(i=12;i<=NF-1;i++) printf \$i"\t"}{print \$NF}' ${covariates}/covariates_HASE_format.txt > pheno_expPC.txt

    mkdir cov_folder
    mkdir pheno_folder

    mv covariate_MDS.txt cov_folder/.
    mv pheno_expPC.txt pheno_folder/.
    """
}

process RunGenRegPcs {

    tag{"Chunk: $Chunk"}

    input:
      path genopath from genotypes_to_genpc
      path covariates from cov_folder
      path phenotypes from pheno_folder
      val studyname from params.studyname
      each Chunk from 1..100

    output:
      file '*_GenRegPcs.txt' into genetic_pcs

    """
    python $baseDir/bin/hase/hase.py \
    -g ${genopath} \
    -study_name ${studyname} \
    -o output \
    -ph pheno_folder \
    -cov cov_folder \
    -th 3 \
    -mode regression \
    -maf 0.01 \
    -node 100 ${Chunk} \
    -cluster "y"

    # calculate degrees of freedom
    N=\$(wc -l pheno_folder/pheno_expPC.txt | awk '{print \$1}')
    N=\$((N-1))
    # nr. of covariates (first column is sample ID but one needs to add SNP here as well)
    # So this is correct (10 gen PCs + SNP = 11)
    N_cov=\$(awk -F' ' '{print NF; exit}' cov_folder/covariate_MDS.txt)
    df=\$((N - N_cov - 1))

    python $baseDir/bin/helperscripts/HaseOutputNumpyAnalyzer.py \
    -i "output/node${Chunk}_*.npy" \
    -df \${df} \
    -o ${Chunk}_GenRegPcs_temp.txt \
    -sref $baseDir/bin/hase/data/1000Gp1v3.ref.gz

    # Filter in only 1e-5
    awk '{if(NR == 1) {print \$0} else {if(\$9 < 1e-5) { print }}}' ${Chunk}_GenRegPcs_temp.txt > ${Chunk}_GenRegPcs.txt
    rm ${Chunk}_GenRegPcs_temp.txt
    """
}

process OrganizeEncodedData {

    input:
      path pd from pd
      path encoded from encoded
      path pd_permuted from pd_permuted
      path encoded_permuted from encoded_permuted
      path mapper from mapper_to_organize
      path genopath from genotypes_to_organise
      val studyname from params.studyname
      file genetic_pcs from genetic_pcs.collectFile(name: 'GenRegPcs.txt', keepHeader: true, sort: true)

    output:
      path '*' into OrganizedFiles

    """
    # empirical
    mkdir -p ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/genotype
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/individuals
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedPhenotypeData
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/empirical/pd_shared

    cp -r ${genopath}/probes ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/
   
    cp ./${encoded}/encode_individuals/* ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/individuals/
    cp ./${encoded}/encode_genotype/* ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedGenotypeData/genotype/
    cp ./${encoded}/encode_phenotype/* ${studyname}_IntermediateFilesEncoded_to_upload/empirical/EncodedPhenotypeData/
    
    cp ./${pd}/*.npy ${studyname}_IntermediateFilesEncoded_to_upload/empirical/pd_shared/

    # permuted
    mkdir -p ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedGenotypeData/genotype
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedGenotypeData/individuals
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedPhenotypeData
    mkdir ${studyname}_IntermediateFilesEncoded_to_upload/permuted/pd_shared

    cp ./${encoded_permuted}/encode_individuals/* ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedGenotypeData/individuals/
    cp ./${encoded_permuted}/encode_genotype/* ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedGenotypeData/genotype/
    cp ./${encoded_permuted}/encode_phenotype/* ${studyname}_IntermediateFilesEncoded_to_upload/permuted/EncodedPhenotypeData/
    
    cp ./${pd_permuted}/*.npy ${studyname}_IntermediateFilesEncoded_to_upload/permuted/pd_shared/

    # Additional files needed for diagnostics
    cp -r ./${mapper} ${studyname}_IntermediateFilesEncoded_to_upload/
    cp -r ${genopath}/SNPQC ${studyname}_IntermediateFilesEncoded_to_upload/
    
    mv GenRegPcs.txt ${studyname}_GenRegPcs.txt
    gzip -f ${studyname}_GenRegPcs.txt
    cp ${studyname}_GenRegPcs.txt.gz ${studyname}_IntermediateFilesEncoded_to_upload/.
    """
}

process ReplaceSampleNames {

    publishDir "${params.outputpath}", mode: 'copy', overwrite: true

    input:
      path OrganizedFiles from OrganizedFiles
      val studyname from params.studyname

    output:
      path "${studyname}_IntermediateFilesEncoded_to_upload" into IntermediateFilesEncodedSampleIdsReplaced_to_upload
      file "*.md5" into md5sumfile

    """
    python $baseDir/bin/helperscripts/replace_sample_names.py -IntFileEnc ${studyname}_IntermediateFilesEncoded_to_upload/empirical/
    python $baseDir/bin/helperscripts/replace_sample_names.py -IntFileEnc ${studyname}_IntermediateFilesEncoded_to_upload/permuted/

    find ${studyname}_IntermediateFilesEncoded_to_upload/ -type f -print0 | xargs -0 md5sum > ${studyname}_OrganizedFiles.md5
    """
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
