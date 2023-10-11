### make gtype file

### LD and Ne
library(vcfR)
library(strataG)
rm(list=ls())
#snp <- read.delim("data/buck.str.geno",sep="\t") 
#snp = snp[,!colnames(snp)%in%c("NA.","CHHS18","CHHS08")]
dat = read.vcfR("data/buckthorn_marineomics_filter.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] = 0
gt2[gt2=="0/1"] = 1
gt2[gt2=="1/1"] = 2
colnames(gt2) = gsub(".bwa.bam","",colnames(gt2))
snp <- t(gt2)# need col = loci; row = ind
### convert into 2 alleles per locus
nloci <- ncol(snp)
############################
### do this code once ######
snp2 <- c()
for (i in 1:nloci)
{
  tmp <- snp[,i]
  key <- data.frame(snp=c(0,1,2),a=c(1,1,2),b=c(1,2,2))
  a <- key$a[match(tmp,key$snp)]
  b <- key$b[match(tmp,key$snp)]
  snp2 <- cbind(snp2,a,b)
}
colnames(snp2) <- sort(c(paste(1:nloci,"a",sep=""),paste(1:nloci,"b",sep="")))
meta = read.csv("data/buckthorn_leaf_meta_META.csv")
meta_sub = meta[match(rownames(snp),meta$ID),]
pop = factor(meta_sub$Site_ID)

dat <- data.frame(rownames(snp),pop,snp2)
colnames(dat) <- c("Ind","Pop",colnames(snp2))
saveRDS(dat,file="data/buck.df")
gi <- df2gtypes(dat, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
saveRDS(gi,file="/data/buck.gtype") ### gtype-formatted data

### performed on larix
gi = readRDS(file="buck.gtype")
# overall
overallHet = heterozygosity(gi, type = "expected") # haplotypic diversity
#### by pop ####
# 1) exp het (haplotypic diversity)
popHet = heterozygosity(gi, type = "expected",by.strata=T) 
# 2) overall PhiST
#popOverall = overallTest(gi,nrep=0)$result["PHIst",1]
# 3) ldNe
ldNe_out0 = ldNe(gi, by.strata=T, maf.threshold = 0,drop.missing=T)
ldNe_out010 = ldNe(gi, by.strata=T, maf.threshold = 0.01,drop.missing=T)
ldNe_out025 = ldNe(gi, by.strata=T, maf.threshold = 0.025,drop.missing=T) 
ldNe_out050 = ldNe(gi, by.strata=T, maf.threshold = 0.05,drop.missing=T) 

# 6) pairwise PhiST
pws <- pairwiseSummary(pairwiseTest(gt,quietly=T,nrep=0))
#phiST_popPair = pws$PHIst; names(phiST_popPair) = paste(pws$strata.1,pws$strata.2,sep="_")

save(overallHet, 
      popHet,
      #popOverall,
      ldNe_out0,
      ldNe_out010,
      ldNe_out025,
      ldNe_out050,
      pws,
      #phiST_popPair, 
      file = "out.RData")
