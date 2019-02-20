phenotype <- as.character(commandArgs(trailingOnly = TRUE))

phenotype <- 30100

path <- "/Users/nbaya/Documents/lab/ukbb-sexdiff/chrX/"

tb <-  read.table(paste0(path,"data/",phenotype,".chrX.tsv.gz"),sep = "\t",header = TRUE)
location <- matrix(unlist(strsplit(as.character(tb$variant),":")), ncol = 4, byrow = TRUE)
  
male <- cbind(location[,c(1:2,4,3)],tb[,c("nCompleteSamples_male","beta_male","se_male","pval_male")])
colnames(male) <- c("CHR","POS","EFFECT_ALLELE","NON_EFFECT_ALLELE","N","BETA","SE","P_VAL")
male$SNP <- tb$variant
male <- male[,c(1,2,9,3:8)]
  
female <- cbind(location[,c(1:2,4,3)],tb[,c("nCompleteSamples_female","beta_female","se_female","pval_female")])
colnames(female) <- c("CHR","POS","EFFECT_ALLELE","NON_EFFECT_ALLELE","N","BETA","SE","P_VAL")
female$SNP <- tb$variant
female <- female[,c(1,2,9,3:8)]

write.table(male,file = paste0(path,"data/",phenotype,".chrX.male.tsv"),
            sep = "\t",row.names = FALSE,quote = FALSE)
write.table(female,file = paste0(path,"data/",phenotype,".chrX.female.tsv"),
            sep = "\t",row.names = FALSE,quote = FALSE)

