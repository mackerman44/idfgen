#' @title Export MSToolKit file
#'
#' @description \code{dumpMSToolKit()} exports a file in MSToolKit format
#'
#' @param x the list of populations to be included in MSToolKit file
#' @param markers the markers to be included in MSToolKit file
#' @param filename what would you like to name the exported MSToolKit file?
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpMSToolKit()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpMSToolKit <- function(x, markers = NULL, filename, replaceBP = TRUE, basePairs = c("A", "C", "G", "T", "-"), replacements = c("1", "2", "3", "4", "5"))
{

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


# Output
  write(paste("", markers, sep = "\t"), filename, ncol = length(markers) * 2, sep = "\t")

  scores <- lapply(pops, twoColScoresAlt, markers)
  output <- do.call(rbind, scores)

  write.table(output, file = filename, append = TRUE, quote = FALSE, sep = "\t ",
            row.names = TRUE, col.names = FALSE)

  cat("A MSToolKit file has been written to", filename, "\n")

}
