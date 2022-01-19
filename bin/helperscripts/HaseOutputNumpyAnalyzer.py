import numpy as np
import timeit
from scipy import stats
import pandas as pd
import os
from collections import OrderedDict
import argparse
import glob
import math
import re

parser = argparse.ArgumentParser(description = "Analyze HASE output Numpy arrays.")

parser.add_argument('-i', '--InputFile', required = True,
                    help = "One or multiple input numpy arrays. Usage of wildcards is supported but then the argument has to be quoted.")

parser.add_argument('-o', '--OutputFile', type = str,
                    required = True,
                    help = "Path to tab-separated output file.")

parser.add_argument('-df', '--DegreesOfFreedom', type = int,
                    required = False,
                    help = "Degrees of freedom for calculating P-values. This should be sample size - number of independent variables in the model (# of covariates + tested) - 1.")

parser.add_argument('-maf', '--MafThreshold', type = float,
                    required = False,
                    help = "MAF threshold for filtering in the SNPs. Has to be 0> and <0.5. Defaults None.", 
                    default = None)

parser.add_argument('-t', '--TThresh', type = float,
                    required = False,
                    help = "Absolute T-value threshold for filtering the results. Defaults None.", 
                    default = None)

parser.add_argument('-phen', '--PhenFilter', required = False, default = None, nargs = '+',
                    help = "One of three: first, file with the list of phenotypes to include. In that case, no header expected. Second, individual phenotype IDs specified and separated by space. Third, specified as None. If not specified, defaults to None.")

parser.add_argument('-snps', '--SnpFilter', required = False, default = None, nargs = '+',
                    help = "One of three: first, file with the list of SNPs/variants to include. In that case, no header expected. Second, individual SNP IDs specified and separated by space. Third, specified as None. If not specified, defaults to None.")

parser.add_argument('-sref', '--SnpRef', required = True,
                    help = "Reference for variant mapper, as specified in HASE documentation. Has to be gzipped and space-delimited.")

# parser.add_argument('-prec', '--RoundingPrecision', type = int,
#                     required = False,
#                     help = "Rounding precision for association statistics. If specified, output file is smaller but statistic precision is reduced. Defaults to None.", 
#                     default = None)

args = parser.parse_args()

# List of numpy arrays:
GeneFilter = args.PhenFilter
SnpFilter = args.SnpFilter
MafTresh = args.MafThreshold
TThresh = args.TThresh
DegreesOfFreedom = args.DegreesOfFreedom
SnpRef = args.SnpRef

if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
    if len(GeneFilter) == 1:
        print GeneFilter
        if re.match(r"[a-zA-Z0-9_.-/]*.txt$", str(GeneFilter).strip('[]"\'')):
            GeneFilter = str(GeneFilter).strip('[]"\'')
            GeneFilter = pd.read_csv(GeneFilter, header = None, delimiter = '\t')
            GeneFilter = set(GeneFilter.iloc[:,0].tolist())
    else:
        GeneFilter = GeneFilter

if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    if len(SnpFilter) == 1:
        print SnpFilter
        if re.match(r"[a-zA-Z0-9_.-/]*.txt$", str(SnpFilter).strip('[]"\'')):
            SnpFilter = str(SnpFilter).strip('[]"\'')
            SnpFilter = pd.read_csv(SnpFilter, header = None, delimiter = '\t')
            SnpFilter = set(SnpFilter.iloc[:,0].tolist())
    else:
        SnpFilter = SnpFilter


# Report some info
print "Input is: " + args.InputFile
print "Output is: " + args.OutputFile
if ((MafTresh is not None) and (MafTresh != 1)):
    print "MAF threshold is: " + str(args.MafThreshold)
if ((TThresh is not None) and (str(TThresh).strip('[]"\'') != 'None')):
    print "T-statistic threshold is: " + str(TThresh)
if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
    print "Output will be filtered to " + str(len(list(GeneFilter))) + " phenotypes."
if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    print "Output will be filtered to " + str(len(list(SnpFilter))) + " SNPs."

input_numpy = glob.glob(args.InputFile)

# Function for P-value calculation
def get_p_value(t_stat, df = None):
    if df is None:
        return stats.norm.sf(np.abs(t_stat))*2
    else:
        return stats.t.sf(np.abs(t_stat), df)*2


# Read in reference
ref = pd.read_csv(SnpRef, compression = 'gzip', sep = ' ')
ref = ref.loc[:, ['ID', 'str_allele1', 'str_allele2']]
if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
    ref = ref[ref['ID'].isin(list(SnpFilter))]

print "SNP reference loaded."

with open(os.path.join(args.OutputFile), 'w') as f:
    f.write('Phenotype\tSNP\tAF\tref_all\talt_all\tt-stat\tbeta\tSE\tP-value\n')
    for filename in input_numpy:
        print "Processing file: " + filename
        batch = np.load(filename).item()
        df = pd.DataFrame(batch)

        # Filtering
        if ((MafTresh is not None) and (MafTresh != 1)):
            print "Filtering in SNPs with MAF> " + str(MafTresh)
            df = df[(df['MAF'] > MafTresh) & (df['MAF'] < (1 - MafTresh))]
        if ((GeneFilter is not None) and (str(GeneFilter).strip('[]"\'') != 'None')):
            print "Gene filter active." 
            df = df[df['phenotype'].isin(list(GeneFilter))]
        if ((SnpFilter is not None) and (str(SnpFilter).strip('[]"\'') != 'None')):
            print "SNP filter active."
            df = df[df['index'].isin(list(SnpFilter))]
        if ((TThresh is not None) and (str(TThresh).strip('[]"\'') != 'None')):
            print "T-statistic filter active."
            df = df[abs(df['t-stat']) > TThresh]

        # Calculate P-values and betas
        Pvals = get_p_value(df["t-stat"], DegreesOfFreedom)
        betas = df['t-stat'] * df['SE']

        # Put into pandas:
        inp = OrderedDict()
        inp["Phenotype"] = df["phenotype"]
        inp["SNP"] = df["index"]
        inp["MAF"] = df["MAF"]
        inp["t-stat"] = df["t-stat"]
        inp["beta"] = betas
        inp["SE"] = df["SE"]
        inp["P-value"] = Pvals

        df_out = pd.DataFrame(inp)

        print df_out[1:5]
        print ref[1:5]

        df_out = pd.merge(df_out, ref, left_on = 'SNP', right_on = 'ID')
        df_out = df_out.rename(columns = {'str_allele1' : 'ref_all', 'str_allele2' : 'alt_all'})

        df_out = df_out[['Phenotype', 'SNP', 'MAF', 'ref_all', 'alt_all', 't-stat', 'beta', 'SE', 'P-value']]

        print df_out[1:5]

        if len(df_out.index) > 0:
            df_out.to_csv(f, sep = "\t", header = False, index = None)
        else:
            print "No rows to include from this array."

print "Output written."
