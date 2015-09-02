# this example is for hugene 1.1 st microarray
# background expression
con2 <- db(pd.hugene.1.1.st.v1)

# Get main probes and neg probes ids. Query below will show which type. 
dbGetQuery(con2, "select * from type_dict")
mainprobes <- dbGetQuery(con2, "select DISTINCT transcript_cluster_id from featureSet where type=1;")
antigm <- dbGetQuery(con2, "select meta_fsetid from core_mps 
inner join featureSet on core_mps.fsetid=featureSet.fsetid where 
featureSet.type='7';") # this can be either 7 or 4, depends on the array
# table(dbGetQuery(con2, "select type from featureSet;")[,1])

# Background expression filter
bkgrdval <- max(apply(exprs(eset)[as.character(antigm[,1]),], 2, quantile,prob=.8))
ind <- genefilter(eset, filterfun(kOverA(280, bkgrdval))) # number of samples 
sum(ind==T)
eset.bkgrdfilter <- eset[ind,] 

# Variance filter
eset.var <- varFilter(eset.bkgrdfilter, var.func=IQR, var.cutoff=0.6, filterByQuantile=TRUE)
