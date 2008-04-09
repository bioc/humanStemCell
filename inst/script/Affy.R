
 library("arrayQualityMetrics")
 library("hgu133plus2.db")

LOAD = FALSE 
if( !LOAD ) {
 hesc = ReadAffy() 
 sN = sampleNames(hesc)
 lcp = lcPrefix(sN)
 lcs = lcSuffix(sN)

 sN = gsub(lcp, "", sN, fixed=TRUE)
 sN = gsub(".CEL", "", sN)

 sampleNames(hesc) = sN

 df = data.frame(SampleID = sN, Diff = rep(c(FALSE, TRUE), c(3,3)))
 rownames(df) = df$SampleID

 vmd = data.frame(labelDescription = c("Sample ID", 
                     "Diff: Differentiated or Not")) 

 phenoData(hesc) = new("AnnotatedDataFrame", 
       data = df, varMetadata = vmd) 

 save(hesc, file="hesc.rda")

} else load("../Affy/hesc.rda")

doQA=FALSE
if(doQA) {
 arrayQualityMetrics(hesc, outdir="QA",
   do.logtransform=TRUE, intgroup="Diff")
}

redoRMA = TRUE
if(redoRMA) {
  hescRMA = rma(hesc)
  fhesc = varFilter(hescRMA)

  ##now we need these genes plus the genes that are present
  load("presGenes.rda")
  allG = union(featureNames(fhesc), presGenes)

  fhesc = hescRMA[allG,]

  fhesc = featureFilter(fhesc)
  
  save(fhesc, file="fhesc.rda")
} else load("../Affy/fhesc.rda")

##fold-change and t-statistics are negative if Undiff has higher expression
##  and positive if Diff has higher expression

  library("limma") 
  design=model.matrix(~fhesc$Diff) 
  hesclim=lmFit(fhesc,design) 
  hesceb=eBayes(hesclim) 

 tab=topTable(hesceb,coef=2,adjust.method="BH", n=6802) 

 tab2 = topTable(hesceb,coef=2,adjust.method="BH", n=15524) 


##do a GSEA/GO analysis here

doGO = FALSE
if(doGO) {
 library(GOstats)

 hgCutoff = 0.01

 ##these correspond to genes that have a lower expression in Diff
 ## (higher in Undiff),  and so get turned off
   
  dIDs = as.character(tab$ID)[tab$t < 0]

  upars = new("GOHyperGParams",
          geneIds= getEG(dIDs, "hgu133plus2"),
          universeGeneIds=getEG(featureNames(fhesc),"hgu133plus2"),
          annotation="hgu133plus2",
          ontology="BP",
          pvalueCutoff=hgCutoff,
          conditional=TRUE,
          testDirection="over")


   uHG =  hyperGTest(upars)
   htmlReport(uHG, file="Upper.html")
   us = summary(uHG)
   us.GO = us$GOBPID
   save(us.GO, file="us.GO.rda")


## these correspond to genes that have a higher expression in Diff, and
## so get turned on

  lIDs = as.character(tab$ID)[tab$t > 0]

 pars = new("GOHyperGParams",
        geneIds= getEG(lIDs, "hgu133plus2"),
        universeGeneIds=getEG(featureNames(fhesc),"hgu133plus2"),
        annotation="hgu133plus2",
        ontology="BP",
        pvalueCutoff=hgCutoff,
        conditional=TRUE,
        testDirection="over")

  lHG =  hyperGTest(pars)
  htmlReport(lHG, file="Lower.html")

  ls = summary(lHG)
  ls.GO = ls$GOBPID
  save(ls.GO, file="ls.GO.rda")
}


