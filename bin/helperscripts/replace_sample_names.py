import pandas as pd
import h5py
import argparse
import numpy as np
import tables
from pandas import HDFStore
import os
from natsort import natsorted

parser = argparse.ArgumentParser(description = "Replace sample IDs in the encoded data to pseudonymized sample IDs. The format is Cohort__SampleNumber. It assumes that data is organised to the folder IntermediateFilesEncoded by previous Nextflow step and it overwrites the files in this folder.")

parser.add_argument('-IntFileEnc', '--IntermediateFilesEncoded', required = True,
                    help = "Folder of encoded data where sample IDs need to be replaced.")

args = parser.parse_args()

print 'Folder with files to share is: ' + args.IntermediateFilesEncoded 

# Get name of the cohort
cohort = os.listdir(args.IntermediateFilesEncoded + '/EncodedGenotypeData/individuals/')
cohort = [x for x in cohort if x.endswith('.h5')]
cohort = cohort[0]
cohort = cohort[:-3]

print 'The name of the cohort is: ' + cohort 

# Encoded genotype individuals file.
individuals = h5py.File(args.IntermediateFilesEncoded + '/EncodedGenotypeData/individuals/' + cohort + '.h5', 'r+')

# Encoded gene expression file
exp = pd.read_csv(args.IntermediateFilesEncoded + '/EncodedPhenotypeData/0_' + cohort + '.csv', sep = "\t", dtype={0:str})
# exp.to_csv(path_or_buf = '/gpfs/hpc/GV/Projects/eQTLGenPhase2/test_HASE_IncompleteData/args.IntermediateFilesEncoded_to_upload2/EncodedPhenotypeData/0_EGCUT_IlluminaArray990.csv', sep = "\t", index = False)
metadata = np.load(args.IntermediateFilesEncoded + '/pd_shared/' + cohort + '_metadata.npy')

# genotye IDs
GenIds = []
for i in range(len(individuals["individuals/table"][:])):
	GenIds.append(individuals["individuals/table"][:][i][1][0])

# phenotype IDs
PhenIds = []
for i in range(exp["id"].shape[0]):
	PhenIds.append(exp["id"][i])

old_ids = list(set(GenIds) | set(PhenIds))
old_ids = natsorted(old_ids)

# Construct new sample IDs.
new_ids = []
for i in range(len(old_ids)):
	new_ids.append(cohort + '__' + str(i))

# Put those into dictionary
old_to_new_dict = dict(zip(old_ids, new_ids))

individuals.close()
## Replace sample IDs in genotype data
# Encoded genotype individuals file.
# Sample ID length is 100 max

assert len(PhenIds) == len(old_to_new_dict), "Sample identifiers are misaligned! expected {} unique samples but observed {}".format(len(PhenIds), len(old_to_new_dict))

# Replace genotype IDs
GenIds = [old_to_new_dict.get(item,item) for item in GenIds]
temp_ind = pd.DataFrame.from_dict({"individual":GenIds})
temp_ind.to_hdf(os.path.join(args.IntermediateFilesEncoded + '/EncodedGenotypeData/individuals/', cohort + '.h5'), key = 'individuals', format = 'table',
				 min_itemsize = 100, complib = 'zlib', complevel = 9)

print('Individuals file adjusted!')

## Replace sample IDs in partial derivatives metadata file
PdIds = list(metadata.item()['id'])
PdIds = [old_to_new_dict.get(item,item) for item in PdIds]

metadata.item()['id'] = np.array(PdIds)
np.save(args.IntermediateFilesEncoded + '/pd_shared/' + cohort + '_metadata.npy', metadata, allow_pickle = True)
print('Sample metadata file for partial derivatives folder is adjusted!')

## Replace sample IDs in encoded expression file
# Replace phenotype IDs
PhenIds = [old_to_new_dict.get(item,item) for item in PhenIds]

for i in os.listdir(args.IntermediateFilesEncoded + '/EncodedPhenotypeData/'):
	exp = pd.read_csv(args.IntermediateFilesEncoded + '/EncodedPhenotypeData/' + i, sep = "\t")
	exp['id'] = PhenIds
	pd.DataFrame(exp).to_csv(args.IntermediateFilesEncoded + '/EncodedPhenotypeData/' + i, sep = "\t", index = False)
	print('Sample names in encoded expression file ' + i + ' are adjusted!')
