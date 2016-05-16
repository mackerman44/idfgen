#' @title Export Structure file
#'
#' @description \code{dumpStructure()} exports a file in Structure format
#'
#' @param x the list of populations to be included in Genepop file
#' @param markers the markers to be included in Genepop file
#' @param filename what would you like to name the exported Genepop file?
#' @param useIndNames would you like to use individual names \code{useIndNames = TRUE} or population name \code{useIndNames = FALSE} as the
#' individual identifiers
#' @param popNumber would you like to include a column in your Structure file that contains the putative population of origin?
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpGenepop()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpStructure <- function(x, markers, filename, useIndNames = TRUE, popNumber = TRUE, replaceBP = TRUE, basePairs = c("A", "C", "G", "T", "-", "0"), replacements = c("1", "2", "3", "4", "5", "-9")) {

# Setup
  if (is.null(markers)) markers <- Markers
# One vs Many Pops
  if(class(x) == "PopList") {
    pops <- lapply(x, get)
  } else {
    pops <- list(x)
  }

# BP Formatting
  if(replaceBP) pops <- lapply(pops, replaceBPs, basePairs, replacements)


# write markers
  write(markers, filename, ncol = length(markers))

# calculate pop column width and create column
  popNames <- unlist(lapply(pops, names))
  indNames <- unlist(lapply(pops, inds))
  indCounts <- unlist(lapply(pops, n))

  if (useIndNames) {
    col1Width <- max(nchar(indNames)) + 1
    column1 <- rep(indNames, rep(2, length(indNames)))
  } else {
    col1Width <- max(nchar(popNames)) + 1
    popColumn <- rep(popNames, indCounts * 2)
  }

  popIndex <- rep(1:length(popNames), indCounts * 2)

# Arrange Scores 2 lines per inds
  structFormat <- function(x, markers) {

    d <- array(dim = c(n(x) * 2, length(markers)))
    d[1:n(x) * 2 - 1, ] <- x$Scores[,markers,1]
    d[1:n(x) * 2, ] <- x$Scores[,markers,2]
    d
  }
  scores <- lapply(pops, structFormat, markers)
  scores <- do.call(rbind, scores)

# Formatting widths of output

  column1 <- format(column1, width = col1Width)
  popIndex <- format(as.character(popIndex), width = 5)
  scores <- format(scores, width = 5)

  if(popNumber) {
    output <- cbind(column1, popIndex, scores)
  } else {
    output <- cbind(column1, scores)
  }

  write.table(output, file=filename, a = TRUE, quote = FALSE, row.names = FALSE, col = FALSE, sep= "")

  cat("A Structure file has been written to", filename, "\n")

}

