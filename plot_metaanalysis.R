phenotype <- as.character(commandArgs(trailingOnly = TRUE))

chrX_path <- "/Users/nbaya/Documents/lab/ukbb-sexdiff/chrX/"

# read meta-analysis results w/ heterogenity and plot
results1 <- read.table(paste0(chrX_path,"data/metaanalysis_results/",phenotype,".metawithhet.1.tsv"),sep = "\t",header = TRUE)

results1$bp <- sapply(strsplit(as.character(results1$MarkerName),":"),function(a) as.numeric(a[2]))

chrmax <- max(results1$bp,na.rm=T)
cols <- c("#004D8B","#269ABC")

results2 <- results1[results1$P.value < 1e-3,]

png(paste0(chrX_path,"plots/","manhattan.",phenotype,".sexdiff.png"),width=183,height=89,res=300,units="mm")
par(mar=c(3.75,3.75,3.2,0.2),mgp=c(2.25,.75,0))

plot(results2$bp,-log10(results2$P.value),
     ylim = c(2.8,max(-log10(results2$P.value)*1.2,-log10(5e-8))),
     col=cols[2],
     pch=20,
     cex=0.6,
     xlab="Chromosome X",
     xaxt="l",
     bty="l",
     ylab="-log10(p)",
     main= paste0("Phenotype: ",phenotype),
     cex.lab=0.9, las=1, cex.axis=0.9)

abline(h=-log10(5e-8),col="red",lty=2,lwd=2)

# end plotting/save
dev.off()


png(paste0(chrX_path,"plots/","manhattan.",phenotype,".hetero.sexdiff.png"),width=183,height=89,res=300,units="mm")
par(mar=c(3.75,3.75,3.2,0.2),mgp=c(2.25,.75,0))

results2 <- results1[results1$HetPVal < 1e-3,]

plot(results2$bp,-log10(results2$HetPVal),
     ylim = c(2.8,max(-log10(results2$HetPVal)*1.2,-log10(5e-8))),
     col=cols[2],
     pch=20,
     cex=0.6,
     xlab="Chromosome X",
     xaxt="l",
     bty="l",
     ylab="-log10(p)",
     main= paste0("Phenotype: ",phenotype, " w/ hetpval"),
     cex.lab=0.9, las=1, cex.axis=0.9)


abline(h=-log10(5e-8),col="red",lty=2,lwd=2)

# end plotting/save
dev.off()
