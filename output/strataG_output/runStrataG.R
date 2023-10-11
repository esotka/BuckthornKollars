
library(strataG)
gi = readRDS(file="buck.gtype")
# overall
overallHet = heterozygosity(gi, type = "expected") # haplotypic diversity
save(overallHet,file="overallHet.Rda")
#### by pop ####
# 1) exp het (haplotypic diversity)
popHet = heterozygosity(gi, type = "expected",by.strata=T) 
save(popHet,file="popHet.Rda")
# 2) overall PhiST
#popOverall = overallTest(gi,nrep=0)$result["PHIst",1]
# 6) pairwise PhiST
pws <- pairwiseSummary(pairwiseTest(gi,quietly=T,nrep=0))
save(pws,file="pws.Rda")
# 3) ldNe
ldNe_out0 = ldNe(gi, by.strata=T, maf.threshold = 0,drop.missing=T)
save(ldNe_out0,file="ldNe_out0.Rda")
#ldNe_out010 = ldNe(gi, by.strata=T, maf.threshold = 0.01,drop.missing=T)
#save(ldNe_out010,"ldNe_out010.Rda")
#ldNe_out025 = ldNe(gi, by.strata=T, maf.threshold = 0.025,drop.missing=T) 
#save(ldNe_out025,"ldNe_out025.Rda")
ldNe_out050 = ldNe(gi, by.strata=T, maf.threshold = 0.05,drop.missing=T) 
save(ldNe_out050,file="ldNe_out050.Rda")

