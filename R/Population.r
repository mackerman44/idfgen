#' @export

#################
# Population.r
#
#
# R S3 CLASS Population
#
# Stucture:
#   x$Scores
#   x$Individuals
#   x$Markers
#   x$IndividualData
#   x$Name
#
#  methods(class="Population")
#


##### R METHODS

  print.Population <- function(x, ...) {
    cat("*** Class Population ***\n")
    cat("Name  :", x$Name, "\n")
    cat("n     :", n(x), "\n")
    cat("m     :", m(x), "\n")
  }
  summary.Population <- function(x, ...) {
    print(x)
    cat("\n*** Sample Scores (Allele 1) ***  - see scores(popName) OR popName$Scores\n\n")
    print(head(x$Scores[,if(m(x) < 6) 1:m(x) else 1:6 ,1]))
    cat("\n*** Sample Individual Data ***  - see metaData(popName) OR popName$IndividualData\n\n")
    print(head(x$IndividualData[,if(length(x$IndividualData) < 4) 1:length(x$IndividualData) else 1:4 ,1]))
  }
  
  c.Population <- function(...) {
    y <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    class(y) <- "PopList"
    y
    
  }

##### Accessors

  n.Population <- function(x) {
    length(x$Individuals)
  }
  
  m.Population <- function(x) {
    length(x$Markers)
  }
  
  inds.Population <- function(x) {
    x$Individuals
  }
  
  markers.Population <- function(x) {
    x$Markers
  }
  
  names.Population <- function(x) {
    x$Name
  }
  
  scores.Population <- function(x) {  # as USER function NOT programming
    oneColScores(x) 
  }

  twoColScores.Population <- function(x, markers = NULL) {
    if (is.null(markers)) markers <- Markers  
    x$Scores[,markers,]
  }
  
  metaData.Population <- function(x) {
    x$IndividualData
  }
##### Custom Functions

  oneColScores.Population <- function(x, markers = NULL) {
     if (is.null(markers)) markers <- x$Markers
     
     d <- paste(x$Scores[,markers,1], x$Scores[,markers,2], sep = "")
     dim(d) <- dim(x$Scores[,markers,, drop=FALSE])[1:2]
     dimnames(d) <- dimnames(x$Scores[,markers,, drop=FALSE])[1:2]
     return(d)
  }
  
  twoColScoresAlt.Population <- function(x, markers = NULL) {
    if (is.null(markers)) markers <- x$Markers
    
    m <- length(markers)
    d <- array(dim = c(n(x), m * 2))
    d[, 1:m * 2 - 1] <- x$Scores[,markers,1]
    d[, 1:m * 2] <- x$Scores[, markers, 2]
    rownames(d) <- inds(x)
    d
    
  
  }
  
  
  
  
  
  
  findDuplicateInds.Population <- function(x, mismatchAllowed, markers = NULL) {
  
    if (is.null(markers)) markers <- x$Markers
    d <- oneColScores(x, markers)
  
  # setup output
    results <- array(NA, dim = c(n(x)*(n(x)-1) / 2, 5))
    colnames(results) <- c("Population", "Ind 1", "Ind 2", "Matches", "Zero Matches")
    results[,1] <- x$Name
    
  # double iteration through all ind pairs
    i <- 0
    for(j in 1:(n(x)-1)){
        for(k in (j+1):n(x)){
            i <- i +1
            results[i,3] <- inds(x)[j]
            results[i,2] <- inds(x)[k] 
            
            matches <- (d[j, markers]==d[k,markers])
            results[i,4] <- sum(matches)
            
            zeroMatches <- d[j,d[j, markers]==d[k,markers]]
            results[i,5] <- length(zeroMatches[zeroMatches=="00"])
        }
    }
  
  # output results   
    out <- results[(as.numeric(results[,4]) >= (length(markers)-mismatchAllowed)), ]
    cat("***", length(out) / 5, "pairs of individuals found in population", x$Name, "with", length(markers) - mismatchAllowed, "or more identical genotypes ***\n") 
    out
                         
  }
  
  
  findNoCalls.Population <- function(x, minNoCall, markers = NULL, ...) {
  
     if (is.null(markers)) markers <- x$Markers
     
     
     d <- oneColScores(x, markers)
      
     truthTable <- d=="00"
     
     zerocount <- apply(truthTable, 1, sum)
     
    
     
     out <- cbind(x$Name, inds(x), zerocount)
     colnames(out) <- c("Population", "Individual", "Zeroes")
     out <- out[as.numeric(out[,3]) >= minNoCall,]
     cat("***", length(out) / 3, "individuals found in population", x$Name, "with", minNoCall, "or more no calls ***\n") 
     out
  
  }
  
  
  removeIndividuals.Population <- function(x, inds, ...) {
        
    rows <- match(inds, inds(x), nomatch = 0)
    rows <- rows[rows != 0]
    if (length(rows) > 0) {
      x$Individuals <- x$Individuals[-rows]
      x$Scores <- x$Scores[-rows, , ]
      x$IndividualData <- x$IndividualData[-rows, ]
      assign(x$Name, x, pos = "IDFGEN_Data")
     
      cat("***", length(rows), "individuals have been removed from population", x$Name, "***\n")
       if(n(x)==0) cat("*** WARNING - All individuals in", x$Name, "have been removed. ***\n")
    
    } 
    
    length(rows)
  }


  findAlleles.Population <- function(x, markers, minorAlleles, ...) {

      a1 <- x$Scores[,markers,1]==rep(minorAlleles, each = n(x))
      a2 <- x$Scores[,markers,2]==rep(minorAlleles, each = n(x))
      counts <- a1 + a2
      
      if(length(markers) > 1) {
        totals <- apply(counts, 1, sum)
      } else {
        totals <- counts
      } 
      indnames <- names(totals[totals>0])
      
      gen <- oneColScores(x)[,markers]
      out <- as.array(cbind(x$Name, inds(x), gen))
      colnames(out) <- c("Population", "Ind", paste(markers, minorAlleles, sep = "_"))
      
      cat("***", length(indnames), "individuals have been found in population", x$Name, "***\n")
      out[indnames,]
      
  }
  
  subPopulation.Population <- function(x, inds, newName, ...) {
  
    newPop <- list()
  
    newPop$Scores <-  x$Scores[inds,,]
    newPop$Individuals <- inds
    newPop$Markers <- x$Markers
    newPop$IndividualData <- x$IndividualData[inds,]
    newPop$Name <- newName
    class(newPop) <- "Population"
    assign(newPop$Name, newPop, pos = "IDFGEN_Data")
    
    a <- newName
    class(a) <- "PopList"
    assign("Populations", c(Populations, a), pos = 1)
    
    cat("*** A new population '", newPop$Name, "' has been created and added to 'Populations' ***\n")
    newPop
    
  
  
  
  }
  
  
  
  
  
  npsScores.Population <- function(x, ...) {
  
    d <- x$Scores
    for (m in markers(x)) {
      d[,m,][d[,m,]==markerAlleles[m,1]] <- 0
      d[,m,][d[,m,]==markerAlleles[m,2]] <- 1
    }
    a <- array(dim = dim(d)[c(1,3,2)])
    a[,1,] <- d[,,1]
    a[,2,] <- d[,,2]
    dimnames(a) <- dimnames(d)[c(1,3,2)]
    a
  }
  
  
  
