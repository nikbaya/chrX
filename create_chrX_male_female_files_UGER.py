#!/broad/software/free/Linux/redhat_6_x86_64/pkgs/anaconda3_2.3.0/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 10:58:58 2018

@author: nbaya
"""

import sys
import pandas as pd
import numpy as np


ph = sys.argv[1]


v3_path = "/psych/genetics_data/projects/ukbb_sexdiff/imputed-v3-results/test/v3"
chrX_path = "/psych/genetics_data/projects/ukbb_sexdiff/imputed-v3-results/test/chrx"


#Males
tb_male = pd.read_csv((v3_path+ph+".imputed_v3.results.male.tsv.gz"), compression='gzip', sep='\t') #read files
chrX_male = tb_male[tb_male.iloc[:]["variant"].str.match('X')][:] #get chrX variants for males
chrX_male = chrX_male.reset_index() #necessary for upcoming concat between chrX_male and a3

a1 = np.asarray(chrX_male.iloc[:,1])
a2 = list(map(lambda variant: str(variant).split(':'), a1))
a3 = pd.DataFrame(np.asarray(a2).reshape((len(a2),4)))

chrX_male2 = pd.concat([a3[[0,1,3,2]],chrX_male], axis = 1).drop(['index','tstat','AC','ytx'], axis =1)
chrX_male2 = chrX_male2.rename(index=str, columns={0: "CHR", 1: "POS", 3: "EFFECT_ALLELE", 2: "NON_EFFECT_ALLELE",
                                     "variant": "SNP", "nCompleteSamples": "N", "beta": "BETA",
                                     "se": "SE", "pval": "P_VAL"})
if len(chrX_male2.columns) > 9:
    chrX_male2 = chrX_male2.drop(['expected_case_minor_AC'],axis=1)
    
chrX_male2.to_csv(chrX_path+ph+".chrX.male.tsv.gz",sep = '\t', compression = 'gzip',index=False)

#Females
tb_female = pd.read_csv((v3_path+ph+".imputed_v3.results.female.tsv.gz"), compression='gzip', sep='\t') #read files
chrX_female = tb_female[tb_female.iloc[:]["variant"].str.match('X')][:] #get chrX variants for females
chrX_female = chrX_female.reset_index() #necessary for upcoming concat between chrX_female and a3

a1 = np.asarray(chrX_female.iloc[:,1])
a2 = list(map(lambda variant: str(variant).split(':'), a1))
a3 = pd.DataFrame(np.asarray(a2).reshape((len(a2),4)))

chrX_female2 = pd.concat([a3[[0,1,3,2]],chrX_female], axis = 1).drop(['index','tstat','AC','ytx'], axis =1)
chrX_female2 = chrX_female2.rename(index=str, columns={0: "CHR", 1: "POS", 3: "EFFECT_ALLELE", 2: "NON_EFFECT_ALLELE",
                                     "variant": "SNP", "nCompleteSamples": "N", "beta": "BETA",
                                     "se": "SE", "pval": "P_VAL"})
if len(chrX_female2.columns) > 9:
    chrX_female2 = chrX_female2.drop(['expected_case_minor_AC'],axis=1)
    
chrX_female2.to_csv(chrX_path+ph+".chrX.female.tsv.gz",sep = '\t', compression = 'gzip',index=False)
