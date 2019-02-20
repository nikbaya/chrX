#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:26:20 2018

@author: nbaya
"""


import os
import glob
import re
import pandas as pd
from subprocess import call
from joblib import Parallel, delayed
import multiprocessing
import sys
import numpy as np


v3_path = "/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/"

#Get saved phenotypes
malefiles = (list(map(os.path.basename,glob.glob(v3_path+"*.male*.gz")))) #restrict to male files to prevent counting phenotype twice
find = re.compile(r"^(.*?)\..*") #regex search term for grabbing all the text before the first period in a string
savedphenotypes = list(map(lambda filename: re.search(find,filename).group(1), malefiles)) #list of all downloaded phenotypes (for me, it gives 78: 77 original samples + 20116_2)

#Get all phenotypes
allphenotypes = pd.Series.tolist(pd.read_table(v3_path+"phenotypes.both_sexes.tsv").iloc[:]["phenotype"]) #list of all phenotypes (male & female)
allphenotypes = pd.DataFrame({'phenotype':allphenotypes})
allphenotypes.to_csv(v3_path+"allphenotypeslist.tsv",sep = "\t")



# TEMPORARY -------------------------------------------------------------------

#savedFiles= (list(map(os.path.basename,glob.glob(chrX_path+"*.gz")))) #restrict to male files to prevent counting phenotype twice
#find = re.compile(r"^(.*?)\..*") #regex search term for grabbing all the text before the first period in a string
#newphenotypes = list(map(lambda filename: re.search(find,filename).group(1), savedFiles)) #list of all downloaded phenotypes (for me, it gives 78: 77 original samples + 20116_2)
#
#nextphenotypes = list(set(savedphenotypes).difference(set(newphenotypes)))
#
#len(nextphenotypes)
# -----------------------------------------------------------------------------


n_cores = multiprocessing.cpu_count()

#old method of extracting chrX
def prev_chrX_from_saved_phenotypes(ph):
    tb_male = pd.read_csv((v3_path+ph+".imputed_v3.results.male.tsv.gz"), compression='gzip', sep='\t') #read files
    tb_female = pd.read_csv((v3_path+ph+".imputed_v3.results.female.tsv.gz"), compression='gzip', sep='\t') 

    chrX_male = tb_male[tb_male.iloc[:]["variant"].str.match('X')][:] #get chrX variants for males
    chrX_female = tb_female[tb_female.iloc[:]["variant"].str.match('X')][:] #get chrX variants for females
    
    
    chrX = pd.merge(chrX_male,chrX_female, on = 'variant',suffixes = ("_male","_female"))
    
    chrX.to_csv(chrX_path+ph+".chrX.tsv.gz",sep = '\t', compression = 'gzip')

#Parallel(n_jobs=n_cores,verbose = 50)(delayed(chrX_from_saved_phenotypes)(ph) for ph in savedphenotypes)

# TEMPORARY -------------------------------------------------------------------
#Parallel(n_jobs=n_cores,verbose = 50)(delayed(chrX_from_saved_phenotypes)(ph) for ph in nextphenotypes)

# -----------------------------------------------------------------------------


#def chrX_from_new_phenotypes(ph):
#    
##        call(["gsutil" ,"cp","gs://ukbb-gwas-imputed-v3-results/export1/"+ph+".**male*",
##              "~/Documents/lab/ukbb-sexdiff/chrX/"])
# 
#
#    call('gsutil ls gs://ukbb-gwas-imputed-v3-results/export1/'+ph+'.**male*', shell=True)
##              "~/Documents/lab/ukbb-sexdiff/chrX/',)
##    call(["paste","<(cat", ph, ".imputed_v3.results.female.tsv.gz","|","zcat", 
##        "|" , "cut -f 1,2,3,5,6,8)", "<(cat", ph,".imputed_v3.results.male.tsv.gz" ,
##        "|", "zcat", "|", "cut", "-f", "1,2,3,5,6,8)", "|", "awk" ,"\'", "NR==1{",
##        "print", "\"variant\",\"n_female\",\"n_male\",\"frq_female\",\"frq_male\",\"beta_female\",\"se_female\",\"p_female\",\"beta_male\",\"se_male\",\"p_male\"",
##        "}NR>1", "&&", "$1==$7{",  "maff=$3/(2*$2);" , "mafm=$9/(2*$8);" , 
##        "if(maff > .05 && maff<.95 && mafm > .05 && mafm < .95){", 
##        "print $1,$2,$8,maff,mafm,$4,$5,$6,$10,$11,$12} }\' | gzip >", ph, ".sexdiff.gz]"])
#
#testph = ['46','47']
#
#for ph in testph:
#    chrX_from_new_phenotypes(ph)

#for ph in set(allphenotypes).difference(set(savedphenotypes)): #for all phenotypes not saved 





# -----------------------------------------------------------------------------
chrX_path = "/Users/nbaya/Documents/lab/ukbb-sexdiff/chrX/data/"

ph = "1757"

#Males
tb_male = pd.read_csv((v3_path+ph+".imputed_v3.results.male.tsv.gz"), compression='gzip', sep='\t') #read files
chrX_male = tb_male[tb_male.iloc[:]["variant"].str.match('X')][:] #get chrX variants for males
chrX_male = chrX_male.reset_index() #necessary for upcoming concat between chrX_male and a3

a1 = np.asarray(chrX_male.iloc[:,0])
a2 = list(map(lambda variant: str(variant).split(':'), a1))
a3 = pd.DataFrame(np.asarray(a2).reshape((len(a2),4)))

chrX_male2 = pd.concat([a3[[0,1,3,2]],chrX_male], axis = 1).drop(['index','tstat','AC','ytx'], axis =1)
chrX_male2.rename(index=str, columns={0: "CHR", 1: "POS", 3: "EFFECT_ALLELE", 2: "NON_EFFECT_ALLELE",
                                     "variant": "SNP", "nCompleteSamples": "N", "beta": "BETA",
                                     "se": "SE", "pval": "P_VAL"})
chrX_male2.to_csv(chrX_path+ph+".chrX.male.tsv.gz",sep = '\t', compression = 'gzip')

#Females
tb_female = pd.read_csv((v3_path+ph+".imputed_v3.results.female.tsv.gz"), compression='gzip', sep='\t') #read files
chrX_female = tb_female[tb_female.iloc[:]["variant"].str.match('X')][:] #get chrX variants for females
chrX_female = chrX_female.reset_index() #necessary for upcoming concat between chrX_female and a3

a1 = np.asarray(chrX_female.iloc[:,0])
a2 = list(map(lambda variant: str(variant).split(':'), a1))
a3 = pd.DataFrame(np.asarray(a2).reshape((len(a2),4)))

chrX_female2 = pd.concat([a3[[0,1,3,2]],chrX_female], axis = 1).drop(['index','tstat','AC','ytx'], axis =1)
chrX_female2.rename(index=str, columns={0: "CHR", 1: "POS", 3: "EFFECT_ALLELE", 2: "NON_EFFECT_ALLELE",
                                     "variant": "SNP", "nCompleteSamples": "N", "beta": "BETA",
                                     "se": "SE", "pval": "P_VAL"})
chrX_female2.to_csv(chrX_path+ph+".chrX.female.tsv.gz",sep = '\t', compression = 'gzip')











