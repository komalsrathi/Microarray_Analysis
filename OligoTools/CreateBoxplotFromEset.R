CreateBoxplotFromEset <- function(eset, listofgenes, var = "treat", nrow = 1, title = "",backgrdline=T)
{
  # check if fData & pData are present
  if(ncol(fData(eset))==0){
    stop("Feature Data absent")
  } 
  
  if(!any(colnames(pData(eset)) %in% var)){
    stop("Phenotypic Data absent or Input variable not recognized")
  }
  
  # get the probeIDs of list of genes & get the corresponding expression
  probeID = as.character(fData(eset)[grep(pattern = paste("^",listofgenes,collapse = "|","$",sep = ""), 
                             x = fData(eset)$Symbol, ignore.case = T),which(colnames(fData(eset)) %in% "ID")])
  if(length(probeID)==1)
  {
    dt = as.data.frame(exprs(eset)[rownames(exprs(eset)) %in% probeID,])
    colnames(dt) = probeID
    dt.melt = melt(t(dt), varnames = c('Var1','Var2'))
  }
  else
  {
    dt = exprs(eset)[rownames(exprs(eset)) %in% probeID,]  
    dt.melt = melt(dt, varnames = c('Var1','Var2'))
  }
  dt.melt = merge(dt.melt, fData(eset), by.x = 'Var1', by.y = 'ID')
  
  # get the variable of interest
  grps = as.character(pData(eset)[,colnames(pData(eset)) %in% var])
  dt.melt = cbind(dt.melt, grps)
  
  # make plot
  p=ggplot(data = dt.melt, aes(x = grps, y = value)) + 
    geom_boxplot(aes(fill = grps), show_guide = F) + 
    #geom_point(aes(color=grps),size=4) + 
    facet_wrap(~Symbol, nrow = nrow) + 
    theme(axis.text.x = element_text(colour = "black", size = 12,angle=90,vjust=0.5),
          axis.text.y = element_text(colour = "black", size = 12), 
          axis.title.x = element_text(colour = "black", size = 12),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(colour = "black", size = 12)) + 
    ggtitle(paste(title,"\n")) + ylab("RMA\n") + xlab(paste('\n',var))
  
  if(backgrdline){
    p=p+geom_hline(aes(yintercept=eset@experimentData@normControls$backgroundExpress), colour="red",linetype="dashed")
  }
  
  return(p)
}
