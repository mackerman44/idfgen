#' @export

####
# PopList.r
#
#
# 

# R FUNCTIONS
print.PopList <- function(x, ...) {
    cat("*** Class PopList ***\n")
    cat("***", length(x), "Populations ***\n")
    print(unclass(x))
    
}

`[.PopList` <- function (x, i)  {
    y <- unclass(x)[i]
    class(y) <- "PopList"
    return (y)
}

c.PopList <- function(...) {
   poplists <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
   objs <- sapply(poplists, get, env = sys.parent())
   allnames <- unlist(objs)
   newPopList <- as.vector(c(allnames))
   class(newPopList) <- "PopList"
   newPopList
}

summary.PopList <- function(x, ...) {
  #print(class(x))
  y <- cbind(n(x))
  colnames(y) <- "Ind Count"
  rownames(y) <- x
  y
}




##### Accessors
scores.PopList <- function(x, f=oneColScores, ... , popCol = FALSE, aslist = FALSE) {
  
  if (popCol) {
      
     pops <- lapply(x, get)
     scores <- lapply(pops, f, ...)
     names <- lapply(pops, names)
     popScores <- mapply(cbind, names, scores, SIMPLIFY=FALSE)
     if(aslist) return(popScores)
     else return(do.call(rbind, popScores))
     
} else {  
  # return scores without pop column
    if(aslist) return(lapply(lapply(x,get),f, ...))
    else return(do.call(rbind, lapply(lapply(x, get), f, ...)))
  
  }
  
}

n.PopList <- function(x, ...) {
  y <- unlist(lapply(lapply(x, get),n))
  names(y) <- x
  y
}
                  
inds.PopList <- function(x, ...) {
  y <- unlist(lapply(lapply(x, get), inds))
  y
}

names.PopList <- function(x, ...) {
  y <- unlist(lapply(lapply(x, get), names))
  y
}


metaData.PopList <- function(x, ...) {

  y <- lapply(lapply(x, get), metaData)
  out <- do.call(rbind, y)
  out
  

}


as.PopList <- function(...) {
  pops <- unlist(list(...))
  class(pops) <- "PopList"
  pops
}



findDuplicateInds.PopList <- function(x, mismatchAllowed, markers = NULL, ...) {

  results <- lapply(lapply(x, get), findDuplicateInds, mismatchAllowed, markers)
  allresults <- t(unlist(sapply(results, t)))
  
  dim(allresults) <- c(5, length(allresults)/5)
  out <- t(allresults)
  cat("***", dim(out)[1], "total individuals found\n")
  colnames(out) <- c("Population", "Ind 1", "Ind 2", "Matches", "Zero Matches")
  rownames(out) <- NULL
  out
}

findNoCalls.PopList <- function(x, minNoCall, markers = NULL) {
  
  results <- lapply(lapply(x, get), findNoCalls, minNoCall, markers)
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  cat("***", dim(out)[1], "total individuals found\n")
  out
}

findAlleles.PopList <- function(x, markers, minorAlleles, ...) {
  results <- lapply(lapply(x, get), findAlleles, markers, minorAlleles)
  allresults <- t(unlist(sapply(results, t)))
  dim(allresults) <- c(length(markers)+2, length(allresults)/(length(markers)+2))
  
   out <- t(allresults)
   out
   cat("***", dim(out)[1], "total individuals found\n")
   colnames(out) <- colnames(results[[1]])
   rownames(out) <- NULL
   out

}

poolPops.PopList <- function(x, newName, ...) {
  
  
  pops <- lapply(x, get)
  d <- lapply(pops, twoColScores)
  meta <- lapply(pops, metaData)

  newPop <- list()
  newPop$Scores <- do.call(abind, list(d, along=1))
  newPop$Individuals <- rownames(newPop$Scores)
  newPop$Markers <- colnames(newPop$Scores)
  newPop$IndividualData <- do.call(rbind, meta) 
  newPop$Name <- newName
  class(newPop) <- "Population"
  assign(newPop$Name, newPop, pos = "IDFGEN_Data")
  
  a <- newName
  print(a)
  class(a) <- "PopList"
  
  assign("Populations", c(Populations, a), pos = 1)
  cat("*** A new population '", newPop$Name, "' has been created and added to 'Populations' ***\n")
  newPop
  

}



