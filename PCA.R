# top 50 variant genes
var <- apply(exprs,1,var)
var <- sort(var,decreasing = T)
var <- data.frame(var = var[1:50])
var <- exprs[rownames(exprs) %in% rownames(var),]
pca <- prcomp(t(var))
scores <- data.frame(Status = pData(eset)$Group, Etiology = pData(eset)$Etiology, Gender = pData(eset)$Gender, 
                    Site = pData(eset)$Site, Sample = pData(eset)$Sample, pca$x)
tmp <- summary(pca)

pc1 <- round(tmp$importance[2,1]*100,2)
pc2 <- round(tmp$importance[2,2]*100,2)

ggplot(data=scores, aes(x=PC1,y=PC2)) + ... # add additional layers
