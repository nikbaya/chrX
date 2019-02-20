

OUTFILE phenotype.metawithhet. .tsv

MARKER   SNP
WEIGHT   N
ALLELE   EFFECT_ALLELE NON_EFFECT_ALLELE
FREQ     EFFECT_ALLELE_FREQ
EFFECT   BETA
STDERR   SE
PVAL     P_VAL

PROCESS /Users/nbaya/Documents/lab/ukbb-sexdiff/chrX/data/phenotype.chrX.male.tsv
PROCESS /Users/nbaya/Documents/lab/ukbb-sexdiff/chrX/data/phenotype.chrX.female.tsv

ANALYZE HETEROGENEITY

QUIT
