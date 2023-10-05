# snmf admix
rm(list=ls())
library(LEA)
dat = read.vcfR("data/buckthorn_marineomics_filter.vcf.gz")
inds = gsub(".bwa.bam","",colnames(dat@gt)[-1])
meta <- read.csv("data/buckthorn_leaf_meta_META.csv")
meta_sub = meta[match(inds,meta$ID),]

### run once
#project = NULL
#project = snmf("data/buck.str.geno", 
#               K=2:20,
#               entropy = TRUE,
#               repetitions = 10,
#               project = "new") 
project = load.snmfProject("data/buck.str.snmfProject")
pdf("output/snmf_cross-val.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
best = which.min(cross.entropy(project, K = 6))
qmatrix <- Q(project,K=6,run=best)
pop_order = factor(meta_sub$Site_ID) # by pop
#pop_order = factor(pop_order,levels=levels(pop_order)[c(6,5,4,2,10,3,1,12,8,11,9,7)])
meta_sub2 = meta_sub[order(pop_order),]
pop_order2 = pop_order[order(pop_order)]
qmatrix_sort = qmatrix[order(pop_order),]

## make a pretty k=6
pdf("output/snmf_k=6.pdf",width=10,height=4)
barplot(t(qmatrix_sort),border=NA,space=0,xlab="Individuals",ylab="Admixture", xlim= c(0,290),ylim=c(-0.18,1.15),
        legend.text= c("1","2","3","4","5","6"),args.legend = list(bty="n",x="right",ncol=1),
        col=c("black","purple","green","orange","red","lightgrey"))
endLine <- as.vector(tapply((1:nrow(qmatrix)),pop_order2,max))
segments(x0=endLine,y0=1,x1=endLine,y1=1.1,col="black",lwd=2)
meanPop <- as.vector(tapply((1:length(pop_order2)),pop_order2,mean))
text(levels(pop_order2),x=meanPop,y=1.07,cex=0.5,srt=90)
dev.off()

kcol <- c("black", #1
          "red",#2
          "darkgreen",#3
          "gainsboro",#4
          "brown",#5
          "deepskyblue",#6
          "yellow",#7
          "dodgerblue4",#8
          "darkorchid",#9
          "burlywood")##10
colorder <- list(
  c(1,2), #ks=2
  c(3,1,2), #ks=3
  c(2,1,4,3), #ks=4
  c(5,4,2,3,1), #ks=5
  c(3,1,6,2,5,4), #ks=6
  c(2,5,3,7,1,6,4), #ks=7
  c(4,7,6,1,2,3,5,8), #ks=8
  c(4,8,3,9,6,5,2,1,7), #ks=9
  c(9,1,2,8,4,6,5,7,10,3)) #ks=10

pdf('output/SNMF.pretty.pdf',width=7,height=4)
#quartz(width=8,height=5)
par(mfrow=c(11,1),mar=c(0,0,0,0),xpd = TRUE)
for (k in 1:9)
{
best = which.min(cross.entropy(project, K = k+1))
qmatrix <- Q(project,K=k+1,run=best)
qmatrix_sort = qmatrix[order(pop_order),]
fig <- barplot(t(qmatrix_sort),col=kcol[colorder[[k]]],space=0,border=NA,xlab="",ylab="",ylim=c(-0.18,1.15),
          names.arg = rep("",nrow(qmatrix_sort)),horiz=F)#,main=substr(ks[k],1,3))
mtext(paste("K",k+1,sep=""),side=2,line=-1.5,at=.5,cex=0.7)
endLine <- as.vector(tapply((1:nrow(qmatrix_sort)),pop_order2,max))
segments(x0=endLine,y0=0,x1=endLine,y1=-.2,col="black",lwd=2)
}
meanPop <- as.vector(tapply((1:length(pop_order2)),pop_order2,mean))
mtext(levels(pop_order2),at=meanPop,side=1,cex=0.5,srt=90)


#xtick<-c(0,3,12,29,35)
#axis(side=1, at=xtick, labels = FALSE)
  
dev.off()


