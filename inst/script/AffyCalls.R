
 library(simpleaffy)

 load("hesc.rda")

 ##this computes present/absent calls using Affy's algorithm

 det = detection.p.val(hesc, tau=0.015, alpha1=0.05, alpha2 = 0.065)

 undiffCalls = det$call[, 1:3]
 diffCalls = det$call[, 4:6]
 ## now select those that were present in at least 2 of the 3 undiff samples
 numPresU = apply(undiffCalls, 1, function(x) sum(x == "P"))

 numPresD = apply(diffCalls, 1, function(x) sum(x=="P"))
 

 presGenes = row.names(undiffCalls)[numPresU >= 2 & numPresD >= 2]


 save(presGenes, file="presGenes.rda")

 presUGenes = row.names(undiffCalls)[numPresU >= 2]
 save(presUGenes, file="presUGenes.rda")

