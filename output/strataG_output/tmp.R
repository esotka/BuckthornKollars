library(strataG)
gi = readRDS(file="buck.gtype")
pws <- pairwiseSummary(pairwiseTest(gi,quietly=T,nrep=1000))
save(pws,file="pws_1000rep.Rda")
