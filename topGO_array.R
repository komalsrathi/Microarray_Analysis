library(topGO)
library(mogene20sttranscriptcluster.db)

# exprset is an Expression Set object
# interesting_genes are the probesets I'd like to perform GO Enrichment on
# first pull out all the probesets

# background data i.e. LIMMA output
bg = read.csv('limma.csv') 

# get background IDs
all_genes <- bg$ID 

# run toGO.func on upregulated genes
godata.upreg <- topGO.func(all_genes = all_genes, interesting_genes = bg[which(bg$fc>0 & bg$adj.P.Val<0.05),2])

# run topGO.func on downregulated genes
godata.dnreg <- topGO.func(all_genes = all_genes, interesting_genes = bg[which(bg$fc<0 & bg$adj.P.Val<0.05),2])

########################################## topGO function ########################################## 
topGO.func <- function(all_genes, interesting_genes)
{
  geneList <- factor(as.integer(all_genes %in% interesting_genes))
  names(geneList) <- all_genes
  
  GOdata <-new ("topGOdata", 
                ontology = "BP", 
                allGenes = geneList, 
                nodeSize = 5, 
                # annot, tells topGO to map from GO terms to "genes"
                annot = annFUN.GO2genes,
                # so annot then calls something to perform this mapping GO2genes, 
                # which is this from the mogene... library
                GO2genes = as.list(mogene20sttranscriptclusterGO2PROBE)
  )
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  allRes <- GenTable(GOdata,
                     classicFisher = resultFisher,
                     classicKS = resultKS,
                     elimKS = resultKS.elim,
                     orderBy = "elimKS", 
                     ranksOf = "classicFisher", 
                     topNodes = 10)
  return(allRes)
}
########################################## topGO function ########################################## 
