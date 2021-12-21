import pandas as pd
import numpy
import os
import argparse
import shutil
import scipy.stats as ss

parser = argparse.ArgumentParser(description = "Convert gene expression or covariate file to HASE format. In case of gene expression matrix it replaces the gene IDs with the ones in the empirical probe mapping file and chunks the table into files <=1,000 phenotypes each. It also enables to inverse normal transform each gene in the input. However, currently batching and inverse normal transformation works only when expression data is fed in, i.e. do not work in case of covariate file is fed in.")

parser.add_argument('-i', '--InputFile', required = True,
                    help = "Path to input gene expression or covariate file. Can be gzipped.")

parser.add_argument('-t', '--Type', required = True,
                    help = "Type of the input. One of exp/cov")

parser.add_argument('-emp', '--EmpiricalProbeMappingFile', required = False,
                    help = "Path to empiricial probe mapping file. Can be gzipped")

parser.add_argument('-gte', '--GenotypeToExpressionFile', type = str,
                    required = False,
                    help = "Path to genotype-to-expression mapping file. Must not contain header.")

parser.add_argument('-o', '--OutputPath', type = str,
                    required = True,
                    help = "Path to output folder. Must not exist.")
                    
parser.add_argument('-n', '--BatchSize', type = int,
                    required = False,
                    help = "Size of the individual batch to divide phenotype table into. HASE seems to fail when batch size is >1000, so defaults to 1000.", 
                    default = 1000)

parser.add_argument('-nc', '--NrOfCovariates', type = int,
                    required = False,
                    help = "Nr. of covariates to include from the first column. This enables to include only first PCs to the covariate file. Defaults to all covariates in the file.")

parser.add_argument('-int', '--InverseNormalTransform', action = 'store_true',
  help = "Boolean specifying whether to apply inverse normal transformation to each gene (or covariate) in the table.")

args = parser.parse_args()

# Messages:
print 'Gene expression or covariate file is:', args.InputFile
print 'Type of the input file:', args.Type
print 'Empirical probe mapping file is:', args.EmpiricalProbeMappingFile
print 'Genotype to expression file is:', args.GenotypeToExpressionFile
print 'Output path is:', args.OutputPath
print 'Nr of phenotypes to add into each batch: ', args.BatchSize
print 'Inverse normal transformation is active: ', args.InverseNormalTransform
if args.Type == "cov":
  if args.NrOfCovariates is None:
    print 'Nr of covariates to include: using all covariates in file'
  else:
    print 'Nr of covariates to include: ', args.NrOfCovariates

# Define function for inverse normal transformation.
def InvNormalTransform(x):
  """
  This function converts numeric input vector to normal distribution by rank-based inverse normal transformation. This function does not accept missing values.
  """
  InvNorm = ss.norm.ppf((ss.rankdata(x) - 0.5) / len(x))
  return InvNorm.tolist()

# Prepare gene expression or covariate matrix
df = pd.read_csv(args.InputFile, sep = "\t")

if args.Type == "exp":
# Prepare emprical probe mapping matrix
  if args.EmpiricalProbeMappingFile != None:
    emp = pd.read_csv(args.EmpiricalProbeMappingFile, sep = "\t")
    emp = emp.drop(axis = 1, labels = "Correlation")

    # Merge
    df["-"] = df["-"].astype('str')
    emp["Probe"] = emp["Probe"].astype('str')
    df2 = pd.merge(df, emp, left_on = "-", right_on = "Probe", how = "inner")
    df2["-"] = df2["Ensembl"]
    df2 = df2.drop(['Ensembl', 'Probe'], axis = 1)
    rownames = list(df2['-'])
    df2.index = rownames
    df2 = df2.drop('-', axis = 1)
  else:
    df2 = df
    df2["-"] = df2["-"].astype('str')
    rownames = list(df2['-'])
    df2.index = rownames
    df2 = df2.drop('-', axis = 1)
  
  # Transpose
  df2 = df2.transpose()

  df2['Pcode'] = df2.index.values
  df2['Pcode'] = df2['Pcode'].astype('str')
  
  # Replace expression codes with genotype codes
  gte = pd.read_csv(args.GenotypeToExpressionFile, sep = "\t", header = None, names = ["Vcode", "Pcode"])
  gte['Pcode'] = gte['Pcode'].astype('str')
  
  df2 = pd.merge(df2, gte, left_on = 'Pcode', right_on = 'Pcode', how = 'inner')
  
  # Rearrange columns
  cols = df2.columns.tolist()
  cols = cols[-1:] + cols[:-2]
  df2 = df2[cols]
  df2 = df2.rename(columns = {"Vcode": "ID"})
  
  print("There are " + str(df2.shape[1]) + " columns in the file.")
  
  # Extract and write out datasets based on specified batch size
  # HASE has bug which seems to cause exeption when batch size is >1000

  cols = df2.columns.tolist()
  
  x = range(1, len(cols))
  
  inp_batches = len(cols) // args.BatchSize + ( len(cols) % args.BatchSize > 0)

  batches = numpy.array_split(numpy.array(x), inp_batches)
  
  # If directory exists, remove it first
  if os.path.isdir(args.OutputPath) == True:
    shutil.rmtree(args.OutputPath)

  os.mkdir(args.OutputPath)
  
  for i in range(0, inp_batches):
    abi_batch = list(batches[i])
    abi_batch.insert(0, 0)
    abi_output = df2[df2.columns[abi_batch]]

    # If inverse normal transformation is required, apply it to each column
    if args.InverseNormalTransform is True:
      inv_normalised = abi_output.drop('ID', axis = 1).apply(InvNormalTransform, axis = 0, result_type = 'broadcast')
      inv_normalised.insert(loc = 0, column = 'ID', value = abi_output['ID'])
      abi_output = inv_normalised
    
    abi_output.to_csv(path_or_buf = args.OutputPath + "/batch" + str(i) + ".txt",  sep = "\t", index = False)
    print(str(i) + " Nr. of traits added: " + str(abi_output.shape[1] - 1))
  
if args.Type == "cov":
  df2 = df
  df2['Pcode'] = df2['-']
  df2['Pcode'] = df2['Pcode'].astype('str')
  
  # Replace expression codes with Vcodes
  gte = pd.read_csv(args.GenotypeToExpressionFile, sep = "\t", header = None, names = ["Vcode", "Pcode"])
  gte['Pcode'] = gte['Pcode'].astype('str')
  
  df2 = pd.merge(df2, gte, left_on = 'Pcode', right_on = 'Pcode', how = 'inner')
  
  # Rearrange columns
  cols = df2.columns.tolist()
  cols = cols[-1:] + cols[:-2]
  cols.pop(1)
  df2 = df2[cols]
  df2 = df2.rename(columns = {"Vcode": "ID"})
  
  if args.NrOfCovariates != None:
    covariates_to_include = args.NrOfCovariates + 1
    df2 = df2.iloc[:, 0:covariates_to_include]
  
  # If directory exists, remove it first
  if os.path.isdir(args.OutputPath) == True:
    shutil.rmtree(args.OutputPath)
    os.mkdir(args.OutputPath)
  df2.to_csv(path_or_buf = args.OutputPath + "/covariates_HASE_format" +".txt",  sep = "\t", index = False)
if args.Type not in ["cov", "exp"]:
  print "Input type has to be on of: exp, cov!"
