#' @title Export GenAlEx file
#'
#' @description \code{dumpGenAlEx()} exports a file in GenAlEx format
#'
#' @param x the list of populations to be included in GenAlEx file
#' @param markers the markers to be included in GenAlEx file
#' @param title what would you like the title of your GenAlEx file to be
#' @param filename what would you like to name the exported GenAlEx file?
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpGenAlEx()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpGenAlEx <- function(x, markers = NULL , title, filename, replaceBP = TRUE, basePairs = c("A", "C", "G", "T", "-"), replacements = c("1", "2", "3", "4", "5")) {

# setup
  if (is.null(markers)) markers <- Markers
  if(class(x) == "PopList") {
    pops <- lapply(x, get)
  } else {
    pops <- list(x)
  }



# Formating
  if(replaceBP) pops <- lapply(pops, replaceBPs, basePairs, replacements)


# Header Rows
  indCounts <- unlist(lapply(pops, n))
  total = sum(indCounts)
  headerRow <- c(length(markers), total, length(pops), indCounts)
  write(headerRow, file = filename, ncolumns = length(headerRow), sep = "\t")

  popNames <- unlist(lapply(pops, names))
  secondRow <- c(title, "", "", popNames)
  write(secondRow, file = filename, ncolumns = length(secondRow), sep = "\t", append = TRUE)

  markerList <- vector()
  markerList[1:length(markers) * 2 - 1] <- markers
  markerList[1:length(markers) * 2] <- ""
  thirdRow <- c("Individual", "Population", markerList)
  write(thirdRow, file = filename, ncolumns = length(thirdRow), sep = "\t", append = TRUE)

# Genetic Data
  outscores <- do.call(rbind, lapply(pops, twoColScoresAlt, markers))
  popcolumn <- rep(popNames, indCounts)

  output <- cbind(popcolumn, outscores)
  write.table(output, file = filename, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


  cat("A GenAlEx file has been written to", filename, "\n")


}
