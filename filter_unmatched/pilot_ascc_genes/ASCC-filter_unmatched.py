# Brendan Reardon
# Van Allen Laboratory
# Dana-Farber Cancer Institute, Harvard Medical School
# 26 February, 2016

# ASCC Paper, filter_unmatched.py

import pandas as pd
import numpy as np
import glob as glob

## ======================
# Import data and subset
## ======================
print 'Running ASCC-filter_unmatched.py'
print 'Importing COSMIC (v75), ExAC, and unmatched samples...'

# We import our dependencies
input_cosmic = '../../storage/reference/cosmic_v75.txt' # This is a frequency table of mutations in COSMIC (v75)
df_cosmic = pd.read_csv(input_cosmic, sep = '\t', low_memory = False)
print '...COSMIC imported.'

# ExAC was prepared as a tab-delimited text file using GATK (https://www.broadinstitute.org/gatk/)
input_exac = '../../storage/reference/exac_GATK.txt'
df_exac = pd.read_csv(input_exac, sep = '\t', low_memory = False)
print '...ExAC imported.'

# And now import our data
allFiles = glob.glob('../../storage/mafs/pilot_snv_unmatched/*.pass.maf')
df_unmatched = pd.DataFrame()
list_ = []

print '...loading mafs...'
for file_ in allFiles:
        df_ = pd.read_csv(file_, index_col = None, sep = '\t', comment = '#', low_memory = False)
        list_.append(df_)
df_unmatched = pd.concat(list_, ignore_index = True)

# Subset for our genes of interest and variant classifications
print '...subsetting to gene list and to variants of interest.'
input_genes = '../../storage/reference/ASCC_genes.txt'
curated_genes = pd.read_csv(input_genes, sep = '\t')
df_unmatched = df_unmatched[df_unmatched['Hugo_Symbol'].isin(curated_genes['gene'])]

variants = ['Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site']
df_variants = pd.DataFrame(variants, columns = ['Variant_Classification'])
df_unmatched = df_unmatched[df_unmatched['Variant_Classification'].isin(df_variants['Variant_Classification'])]

# We require 14 reads covering a site for declaring a site to be adequately covered for mutation calling
df_unmatched['read_depth'] = df_unmatched['t_alt_count'] + df_unmatched['t_ref_count']
df_unmatched = df_unmatched[df_unmatched['read_depth'] >= 14]

df_unmatched.index = range(0, len(df_unmatched)) # Reset the index of the dataframe

## ======================
# Annotate with Cosmic
## ======================
print 'Annotating', str(len(df_unmatched)), 'mutations with COSMIC(v75)...'

# Add a blank Cosmic(v75) counts column to dataframe. In this column, we will record how often this mutation is observed in the COSMIC database.
df_unmatched['Cosmic(v75) Counts'] = 0
for i in range(0, len(df_unmatched)):
    if len(df_cosmic[(df_cosmic['Gene name'] == df_unmatched.ix[i,'Hugo_Symbol']) & (df_cosmic['Mutation AA'] == df_unmatched.ix[i,'Protein_Change'])]) == 0:
        df_unmatched.ix[i,'Cosmic(v75) Counts'] = np.nan
    else:
        df_unmatched.ix[i,'Cosmic(v75) Counts'] = df_cosmic[(df_cosmic['Gene name'] == df_unmatched.ix[i,'Hugo_Symbol']) & (df_cosmic['Mutation AA'] == df_unmatched.ix[i,'Protein_Change'])]['count'].get_values()[0]
    # We print a small progress report for the user
    if i == round(len(df_unmatched)/4) or i == round(len(df_unmatched)/2) or i == round(len(df_unmatched)/1.31) or i == (len(df_unmatched) - 1):
        print('...' + str(i) + ' of ' + str(len(df_unmatched)) + ' iterations complete. \n')
df_unmatched['Cosmic(v75) Counts'] = df_unmatched['Cosmic(v75) Counts'].fillna(0)

print '...Cosmic annotation complete.'

## ======================
# Filter with ExAC
## ======================
print 'Filtering with ExAC....'

df_tmp = df_unmatched.copy(deep =True)
df_tmp_indexes = np.transpose(df_tmp.index.values.tolist())

df_germ = pd.DataFrame([], columns = df_unmatched.columns.values.tolist())
for i in range(0, len(df_tmp)):
    ind = df_tmp_indexes[i]
    tmp_alt = df_tmp.loc[ind, 'Tumor_Seq_Allele2']
    tmp_row = df_tmp.loc[[ind]]
    tmp_exac = df_exac[(df_exac['CHROM'] == df_tmp.loc[ind,'Chromosome']) & (df_exac['POS'] == df_tmp.loc[ind,'Start_position'])]
    if len(tmp_exac) != 0:
        tmp_exac_index = tmp_exac.index.values.tolist()[0]
        if (tmp_exac['ALT'].str.contains(tmp_alt)[tmp_exac_index] == True):
            df_germ = df_germ.append(tmp_row)
        # We print a small progress report for the user
    if i == round(len(df_tmp)/4) or i == round(len(df_tmp)/2) or i == round(len(df_tmp)/1.31) or i == (len(df_tmp) - 1):
        print('...' + str(i) + ' of ' + str(len(df_tmp)) + ' iterations complete. \n')

df_exac_passed = df_unmatched.drop(df_germ.index)
df_exac_failed = df_germ

print('...ExAC filter complete.')

print('Exporting results from filtering with ExAC...')
df_exac_passed.to_csv('df_exac_passed.txt', sep = '\t', index = False)
df_exac_failed.to_csv('df_exac_failed.txt', sep = '\t', index = False)

## ======================
# Recover with COSMIC
## ======================
print 'Readding variants that appear in Cosmic at least 3 times...'

# We subset the failed mutations for those that had at least 3 occurences in Cosmic
df_recover = df_exac_failed[df_exac_failed['Cosmic(v75) Counts'] >= 3]

# and concatenate them together with the mutations that passed
df_filter_pass = pd.concat([df_exac_passed, df_recover], ignore_index = True)
df_filter_fail = df_unmatched.drop(df_filter_pass.index)

# And Export!
print 'Exporting overall results...'
df_filter_pass.to_csv('df_filter_pass.txt', sep = '\t', index = False) # This is the primary file that we will be using as we move forward
df_filter_fail.to_csv('df_filter_failed.txt', sep = '\t', index = False) 

print 'filter_unmatched.py complete!'
