# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:34:26 2024

@author: kcoykendall
"""

# -*- coding: utf-8 -*-
# kcoykendall July 2023
##Counts alleles of the sex marker
##Sums reads
##Counts forward primer instances
##Calculates Allele1/Allele2 ratios
## This case is a ZW marker. Allele1=T8, Allele2=T9
## Females typically would be T8/T9 and males T9/T9
## Pay attention to the type of fastq you are working with (initial, re-run, fill-in). The
###sample names will have different prefixes you'll want to remove
#Requires a pheno file of sampleID in column 1 and pheno sex in column 2

import os
import pandas as pd
import numpy as np

os.getcwd()
os.chdir('./Desktop/Test/')
###variables to consider
directory = r'./' 
Allele1="GGTTTTTTTTC"
Allele2="GGTTTTTTTTTC"
#forward sequence
forseq = 'GCTGTCAGATAAATGTAGTGAAACAAA'
#Tab delimited file, with header, with sample name in first column, phenotypic sex in second column
pheno_file = "ZW_pheno_L0450.txt"
#name of file that will be exported
ExportZWinfo='ZW_info_L0450.csv'

### counting forward primer seq
forward = []
for filename in os.listdir(directory):
  if filename.endswith(".fastq"):
    f=open(os.path.join(directory, filename),'rt')
    num_for = 0
    for line in f:
        if forseq in line:
            num_for += 1
    forward.append((filename, num_for))
    print('finished with', filename, 'forward primer counts')  

for_data = pd.DataFrame(forward, columns=['Sample', 'forward_primer'])
for_data['forward_primer'] = for_data['forward_primer'].astype(float)

###counting reads with forward sequence AND allele1  
Allele1List = []
for filename in os.listdir(directory):
  if filename.endswith(".fastq"):
    f=open(os.path.join(directory, filename),'rt')
    Allele1Count = 0
    for line in f:
        if forseq in line and Allele1 in line:
            Allele1Count += 1
    Allele1List.append((filename, Allele1Count))
    print('finished with', filename, 'Allele1 counts')    
Allele1_data = pd.DataFrame(Allele1List, columns=['Sample', 'Allele1'])
Allele1_data['Allele1'] = Allele1_data['Allele1'].astype(float)

###counting reads with forward sequence AND allele2 
Allele2List = []
for filename in os.listdir(directory):
  if filename.endswith(".fastq"):
    f=open(os.path.join(directory, filename),'rt')
    Allele2Count = 0
    for line in f:
        if forseq in line and Allele2 in line:
            Allele2Count += 1
    Allele2List.append((filename, Allele2Count))
    print('finished with', filename, 'Allele2 counts')

Allele2_data = pd.DataFrame(Allele2List, columns=['Sample', 'Allele2'])
Allele2_data['Allele2'] = Allele2_data['Allele2'].astype(float)

###put results together in one dataframe, add sum and ratio
BothAll=Allele1_data.merge(Allele2_data, left_on='Sample', right_on='Sample')
ZW_counts=BothAll.merge(for_data, left_on='Sample', right_on='Sample')

###change 0's to 0.00001 so there won't be any division by 0
ZW_counts=ZW_counts.replace(0.0, 0.00001)
ZW_counts['Sum'] = ZW_counts['Sum'] = ZW_counts.loc[:, 'Allele1':'Allele2'].sum(1)
ZW_counts['Ratio'] = ZW_counts['Allele1']/ZW_counts['Allele2']
# Note that in the case where there are zero reads for both alleles, their allele
#counts were changed to 0.00001 so the value of Allele1/Allele2 will be 1. Next line changes
#value in Allele1/Allele2 back to 0
ZW_counts['Ratio'] = np.where((ZW_counts['Allele1'] == 0.00001) & (ZW_counts['Allele2'] == 0.00001), '0.0', ZW_counts['Ratio']) 
ZW_counts['Ratio'] = ZW_counts['Ratio'].astype(float)

###Genotype the ZW marker based on set of conditions similar to our 
###usual genotyping thresholds EXCEPT in this case, if the ratio is over 0.75
##that will be called a female
###if the sum of allele1 and allele2 counts is >10, otherwise Results=TooLow
###if the ratio of Allele1 to Allele 2 is >0.75, calls Female
###if the ratio of Allele1 to Allele 2 is <0.2, calls Male
###All other values for Allele1/Allele2 (0.21 - 0.74) called 'WeirdRatio' and not genotyped 
###NOTE! My conditions allow for females to have only Allele1 and no allele2 (Allele1/Allele1)

conditions = [
    ZW_counts['Sum'] <= 10.0, #no genotype if sum of counts are <=10
    ZW_counts['Ratio'] >= 10, #if the ratio of Allele1 to Allele2 is 10 and above, WW, Female
    (ZW_counts['Ratio'] <= 3.0) & (ZW_counts['Ratio'] >=0.70), #if the ratio of Allele1 to Allele2 is between 0.7 and 3.0, Female
    ZW_counts['Ratio'] <= 0.2]  #if ratio of Allele1 to Allele2 is less than 0.2, Male
#Note that any 8T/9T values between 0.2 and 0.74 will be scored as 'NG' (default value below)
choices = ['TooLow', 'WW', 'ZW', 'ZZ']
ZW_counts['Results'] = np.select(conditions, choices, default='WeirdRatio')

ZW_counts['GenoSex'] = ""
ZW_counts['GenoSex'] = np.where((ZW_counts['Results'] == 'ZW'), 'Female', ZW_counts['GenoSex'])
ZW_counts['GenoSex'] = np.where((ZW_counts['Results'] == 'ZZ'), 'Male', ZW_counts['GenoSex'])
ZW_counts['GenoSex'] = np.where((ZW_counts['Results'] == 'WW'), 'Female', ZW_counts['GenoSex'])
ZW_counts['GenoSex'] = np.where((ZW_counts['Results'] == 'TooLow'), 'Unknown', ZW_counts['GenoSex']) 
ZW_counts['GenoSex'] = np.where((ZW_counts['Results'] == 'WeirdRatio'), 'Unknown', ZW_counts['GenoSex'])

##Eagle Fish Genetics Lab affixes 1 of 4 prefixes to fastqs when genotyping: initial, r1, f1, qc
##deleting the prefix at the beginning of the sample names and ".fastq" at end
ZW_counts['Sample'] = ZW_counts['Sample'].str.replace('.fastq', '')
ZW_counts['Sample'] = ZW_counts['Sample'].str.replace('f1', '')
ZW_counts['Sample'] = ZW_counts['Sample'].str.replace('r1', '')
ZW_counts['Sample'] = ZW_counts['Sample'].str.replace('qc', '')
ZW_counts['Sample'] = ZW_counts['Sample'].str.replace('initial', '')

## get the phenotypic sexes. 
pheno_sexes = pd.read_csv(pheno_file, sep='\t')
ZW_info = pd.merge(ZW_counts, pheno_sexes, on="Sample")
ZW_info['Concordance'] = np.where((ZW_info['GenoSex']==ZW_info['Pheno']), 'Concordant', 'Discordant')
ZW_info['Concordance'] = np.where((ZW_info['GenoSex']=="Unknown"), 'NA', ZW_info['Concordance'])
ZW_info['Concordance'] = np.where((ZW_info['Pheno']=="Unknown"), 'NA', ZW_info['Concordance'])

## export table
ZW_info.to_csv(ExportZWinfo, sep=',', na_rep='', float_format=None, columns=None, header=True)

