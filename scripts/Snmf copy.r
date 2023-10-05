### make snmf
library(vcfR)
library(poppr)
library(hierfstat)
library(LEA)
rm(list=ls())

setwd("C:/Users/nmkol/My Drive/Hughes_Lab/Kollars_HughesLab_Postdoc_Work/Buckthorn/Nicole_analysis")

dat = read.vcfR("buckthorn_marineomics_filter.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] = 0
gt2[gt2=="0/1"] = 1
gt2[gt2=="1/1"] = 2
#table(is.na(gt2))
#  FALSE    TRUE 
#1960713   18607 
gt2_genoYN = is.na(gt2)
hist(1-rowSums(gt2_genoYN)/ncol(gt2)) # loci with missing data
table(1-rowSums(gt2_genoYN)/ncol(gt2)) #==> 2315 SNPs with no NAs
hist(1-colSums(gt2_genoYN)/nrow(gt2)) # ind with missing data
table(1-colSums(gt2_genoYN)/nrow(gt2)) #==> 42 ind with no NAs


snp.hf <- df2genind(t(gt2),ncode=1)
snp.hf2 <- genind2hierfstat(dat=snp.hf,pop = sample.meta.order$Site_ID)
write.struct(snp.hf2,ilab=sample.meta.order$Sample.geno,pop=sample.meta.order$Site_ID,fname="buck.str",MARKERNAMES = F,MISSING=9)
### remove NAs in first column - BUT WHERE DO THE NAS COME FROM???

struct2geno("buck.str",ploidy=2,FORMAT = 2) 

project = NULL
project = snmf("buck.str.geno", 
               K=2:10,
               entropy = TRUE,
               repetitions = 1, # need to include more
               project = "new") #19 is the max populations sampled
project = load.snmfProject("buck.str.snmfProject")
pdf("snmf_cross-val.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

best = which.min(cross.entropy(project, K = 3))
qmatrix <- Q(project,K=3,run=best)
pop = substr(colnames(gt2),1,2)
pop_order = factor(pop) # by pop

qmatrix_sort = qmatrix[order(pop_order),]
pop_order2 = pop_order[order(pop_order)]

pdf("snmf.pdf",width=10,height=4)
barplot(t(qmatrix_sort),border=NA,space=0,xlab="Individuals",ylab="Admixture")#, xlim= c(0,235),
        #legend=TRUE,legend.text= c("1","2","3","4","5","6"),args.legend = list(bty="n",x="right",ncol=1),
        #col=c("black","purple","green","orange","red","lightgrey"))
endLine <- as.vector(tapply((1:length(pop_order2)),pop_order2,max))
segments(x0=endLine,y0=-1,x1=endLine,y1=1,col="red",lwd=2)
meanPop <- as.vector(tapply((1:length(pop_order2)),pop_order2,mean))
mtext(levels(pop_order2),at=meanPop,cex=0.5)
dev.off()



