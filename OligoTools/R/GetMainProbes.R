#'Annotate a ExpressioSet lifted from Affycoretools
#' @param eset An ExpressionSet
#'
#'
#'

GetMainProbes <- function(input){
  {
    if (is(input, "ExpressionSet")) {
      pdinfo <- annotation(input)
      if (length(grep("^pd", pdinfo)) != 1) 
        stop(paste("The file", pdinfo, "does not appear to have been processed using", 
                   "the oligo package.\nIn this case the argument to this function should", 
                   "be the name of the correct pd.info package (e.g., pd.hugene.1.0.st.v1.\n"), 
             call. = FALSE)
    }
    else {
      if (is.character(input)) 
        pdinfo <- input
      else if (!is.character(input)) 
        stop(paste("The input argument for this function should either be an ExpressionSet", 
                   "that was generated using oligo, or the name of the pd.info package", 
                   "that corresponds to your data.\n"), call. = FALSE)
    }
    require(pdinfo, character.only = TRUE)
    con <- db(get(pdinfo))
    types <- dbGetQuery(con, paste("select distinct meta_fsetid, type from featureSet inner join core_mps", 
                                   "using(fsetid);"))
    dbDisconnect(con)
    if (is(input, "ExpressionSet")) {
      types <- types[match(featureNames(input), types[, 1]), 
                     ]
      ind <- types[, 2] %in% 1
      return(input[ind, ])
    }
    else return(types)
  }
  

}
