### AMOVA on adults
### method 1 = pegas::vegan on allele counts
library(vcfR)
rm(list=ls())
## Genotype calls
meta = read.csv("data/buckthorn_leaf_meta_META.csv")
vcf = read.vcfR("data/buckthorn_marineomics_filter.vcf.gz")
gt <- extract.gt(vcf,return.alleles=F)

gt[gt=="0/0"] = 0
gt[gt=="0/1"] = 1
gt[gt=="1/1"] = 2
colnames(gt) = gsub(".bwa.bam","",colnames(gt))
gt2 = as.data.frame(gt)
meta_sub = meta[match(colnames(gt2),meta$ID),]

## AMOVA

st.d <- dist(t(gt2))
groups = factor(meta_sub$MetaPop)
pop = factor(meta_sub$Site_ID)
m1 <- amova(st.d~groups/pop)
sink("output/amova_out.txt")
print(m1)
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))
sink()


