#' @title Export Genepop file
#'
#' @description \code{dumpGenepop()} exports a file in Genepop format
#'
#' @param x the list of populations to be included in Genepop file
#' @param markers the markers to be included in Genepop file
#' @param title what would you like the title of your Genepop file to be
#' @param filename what would you like to name the exported Genepop file?
#' @param popNames would you like to use population names \code{popNames = TRUE} or individual names \code{popNames = FALSE} as individual
#' identifiers in the Genepop file
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

dumpGenepop <- function(x,
                        markers = NULL,
                        title,
                        filename,
                        popNames = FALSE,
                        replaceBP = TRUE,
                        basePairs = c("A", "C", "G", "T", "-"),
                        replacements = c("1", "2", "3", "4", "5")) {


# Markers
  if (is.null(markers)) markers <- Markers
  snps <- markers[markers %in% MarkersSNP]
  usats <- markers[markers %in% MarkersUSAT]


# One vs Many Pops
  if(class(x) == "PopList") {
    pops <- lapply(x, get)
  } else {
    pops <- list(x)
  }

# Formating
  if(replaceBP) pops <- lapply(pops, replaceBPs, basePairs, replacements)

  formatScores <- function(x, snps, usats, popNames) {

    x$Scores[,snps,] <- paste("0", x$Scores[,snps,], sep = "")
    x$Scores[,usats,][x$Scores[,usats,]==0] <- "000"
    gen <- oneColScores(x)
    if (popNames==FALSE) {
      out <- cbind(paste(inds(x), ",", sep = ""), gen[,markers])
    } else {
      out <- cbind(paste(names(x), ",", sep = ""), gen[,markers])
    }
    out
  }

  genepopScores <- lapply(pops, formatScores, snps, usats, popNames)


# Output
  write(c(title, markers), file = filename, ncolumns = 1, sep = "")

  writePop <- function(x) {
    write("pop", file = filename, ncolumns = 1, sep = "", append = TRUE)
    write.table(x, file = filename, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
  }

  lapply(genepopScores, writePop)
  cat("***", filename, "has been written. ***\n")


}
