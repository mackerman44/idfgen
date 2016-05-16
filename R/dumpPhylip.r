#' @title Export Phylip file
#'
#' @description \code{dumpPhylip()} exports a file in Phylip format
#'
#' @param x the list of populations to be included in Phylip file
#' @param markers the markers to be included in Phylip file
#' @param filename what would you like to name the exported Phylip file?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpPhylip <- function(x, markers = NULL, filename) {
#### needs code refactoring


# Setup
  if (is.null(markers)) markers <- Markers
# One vs Many Pops
  if(class(x) == "PopList") {
    pops <- lapply(x, get)
  } else {
    pops <- list(x)
  }

# write header
  write(c(length(x), length(markers)), filename, ncol = 2)



# get all scores for all populations to build allele table
  structFormat <- function(x, markers) {

    d <- array(dim = c(n(x) * 2, length(markers)))
    d[1:n(x) * 2 - 1, ] <- x$Scores[,markers,1]
    d[1:n(x) * 2, ] <- x$Scores[,markers,2]
    d
  }
  scores <- lapply(pops, structFormat, markers)
  allScores <- do.call(rbind, scores)
  colnames(allScores) <- markers
# build allele table

  a <<- allScores
  alleleTable <- lapply(alply(allScores, 2, unique), sort)
  removeZeros <- function(x) {
    x[x!=0]
  }
  alleleTable <- lapply(alleleTable, removeZeros)

# write line 2
  numOfAlleles <- sapply(alleleTable, length)
  write(numOfAlleles, filename, ncol = length(numOfAlleles), sep = " ", append = TRUE)

# calculate freqs for each pop
  popNames <- unlist(lapply(pops, names))
  for (i in 1:length(popNames)) {

    write(format(popNames[i], width = 10), filename, ncol = 1, append = TRUE)
    scores <- rbind(pops[[i]]$Scores[, markers, 1], pops[[i]]$Scores[, markers, 2])
    for (j in 1:length(alleleTable)) {
      alleleCounts <- ""

      for (a in 1:length(alleleTable[[j]])) {
         alleleCounts[a] <- length(grep(alleleTable[[j]][a], scores[,j]))
      }
      freqs <-  as.numeric(alleleCounts) / sum(as.numeric(alleleCounts))
      write(sprintf("%1.4f", freqs), filename, ncol = length(freqs), sep = " ", append = TRUE)
    }
  }
  cat("A Phylip file has been written to", filename, "\n")
  cat("\n** currently dumpPhylip 'infile' output does not work in program Phylip. Further testing is needed to resolve. **\n")

}

