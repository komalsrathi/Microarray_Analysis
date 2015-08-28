# basic libraries to load
library(oligo)
library(reshape)
library(ggplot2)
library(gridExtra)
library(limma)

# libraries based on the genechip
library(mogene20sttranscriptcluster.db)

# set working directory

########### load sample/phenotype data ################
s <- read.csv("/data/samples.csv") #usually contains the id, type/treat of sample
grps <- s$type  # type of samples
names <- s$id   # id of samples
########### load sample/phenotype data ################


################## get normalized counts #########################
(celFiles <- list.celfiles(getwd())) # get list of cel files
affyGeneFS <- read.celfiles(celFiles,verbose=TRUE) # read celFiles
con2 <- db(pd.mogene.2.0.st) # db connection
eset <- rma(affyGeneFS, target='core') # # rma normalization GeneFeatureSet
eset <- rma(affyGeneFS) # # rma normalization ExpressionFeatureSet

expr.unnorm <- exprs(affyGeneFS) # contains unnormalized counts of all genes
expr <- exprs(eset) #contains normalized counts of all genes
################## get normalized counts #########################



########################## boxplots ############################
# boxplot of raw un-normalized counts
melted <- melt(t(expr.unnorm))
melted <- cbind(melted,grps,names)
ggplot(data=melted, aes(x=as.character(names),y=log2(value),colour=grps))+ geom_boxplot() + theme(axis.text.x =element_text(angle = 90, vjust = 0.5,size=14,colour="black"),axis.text.y=element_text(size=14,colour="black"), axis.title.x=element_blank(),axis.title.y=element_text(size=14),title=element_text(size=14,face="bold")) + ggtitle("Boxplot:\nRaw log2 intensities")

# boxplot of normalized counts
melted <- melt(t(expr))
melted <- cbind(melted,grps,names)
ggplot(data=melted, aes(x=as.character(names),y=value,colour=grps))+ geom_boxplot() + theme(axis.text.x =element_text(angle = 90, vjust = 0.5,size=14,colour="black"),axis.text.y=element_text(size=14,colour="black"), axis.title.x=element_blank(),axis.title.y=element_text(size=14),title=element_text(size=14,face="bold")) + ggtitle("Boxplot:\nRMA Normalized intensities")

grid.arrange(p1,p2,ncol=2)
########################## boxplots ############################



##################### PCA ######################################
var <- apply(t(expr),1,var)
var <- sort(var,decreasing=T)
var <- as.data.frame(var[1:500]) # top 500 variant genes

top.var <- dat[,which(colnames(expr) %in% row.names(var))] # get expression for the same

pca <- prcomp(top.var) 
scores <- data.frame(pca$x[,c(1:2)]) # get PC1 & PC2
scores$grps <- grps
scores$names <- names
ggplot(data=scores, aes(x=PC1,y=PC2)) + geom_text(aes(label=names,color=factor(grps)),cex=5) + xlab("PC1") + ylab("PC2") + ggtitle("PCA Top 500 Variant Genes\n")
##################### PCA ######################################



##################### Filter out unexpressed genes #############
# Get main probes and neg probes ids. Query below will show which type. 
dbGetQuery(con2, "select * from type_dict")
mainprobes <- dbGetQuery(con2, "select DISTINCT transcript_cluster_id from featureSet where type=1;") # main
negprobes  <- dbGetQuery(con2, "select DISTINCT transcript_cluster_id from featureSet where type=4;") # antigenomic
negeset <- expr[which(rownames(expr) %in% negprobes$transcript_cluster_id),] 
negavg <- mean(negeset) 

expr.filtered <- expr[apply(expr,1,FUN=function(x) all(x > negavg)),] # all values > than negavg are kept
colnames(expr.filtered) <- s$id
expr.filtered # contains normalized counts of expressed genes
##################### Filter out unexpressed genes #############



################################### DEGs #########################
# differential gene expression using LIMMA
design <- model.matrix(~0 + type,data=s) #change depending on sample info
colnames(design) <- levels(s$type)
fit <- lmFit(expr.filtered, design)
cont.matrix <- makeContrasts(grp=mut-ctrl,
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2,coef=1,number=Inf,p.value=0.05,adjust.method="BH") # DEGs

# calc fold change from log2FC
results$FC<-2^(abs(results$logFC))
results$FC<-ifelse(results$logFC<0,results$FC*(-1),results$FC*1)
################################### DEGs #########################



################################## Annotate DEGs #################
# get geneNames & ensemblIds for the genechip
Annot <- data.frame(SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=","), 
                    ENSEMBLID=sapply(contents(mogene20sttranscriptclusterENSEMBL), paste, collapse=","))

# add annotation to results
results.ann <- merge(results,Annot,by="row.names") 
################################## Annotate DEGs #################

