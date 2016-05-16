#' @title Export GSI baseline file
#'
#' @description \code{dumpBaseline()} exports a gsi baseline file in BAYES, SPAM, or gsi_sim format
#'
#' @param x the list of populations to be included in gsi baseline file
#' @param markers the markers to be included in gsi baseline file
#' @param fileType the format of the exported gsi baseline file. Can be set to "SPAM", "BAYES", or "gsi_sim"
#' @param filename what would you like to name the exported gsi baseline file
#' @param findMarkers i don't remember what this does, just leave as the default \code{NULL}
#' @param replaceWith something related to \code{findMarkers}. Also leave as the default \code{NULL}
#' @param removePrefix this can be used to trim the marker names for "SPAM"
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpBaseline()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpBaseline <- function(x, markers, fileType, filename, findMarkers = NULL, replaceWith = NULL, removePrefix = NULL, replaceBP = TRUE, basePairs = c("A", "C", "G", "T", "-"), replacements = c("1", "2", "3", "4", "5")) {
###############################################################################
#
#   Outputs baseline in BAYES, SPAM, or gsi_sim format
#
#       updated 12/29 to add gsi_sim by EG
#       updated 2/21 for version 2.1 by EG

# setup
  if (is.null(markers)) markers <- Markers
  if(class(x) == "PopList") {
    pops <- lapply(x, get)
  } else {
    pops <- list(x)
  }
  m <- length(markers)

# Formating



  if (fileType == "gsi_sim") {

    if(replaceBP) pops <- lapply(pops, replaceBPs, basePairs, replacements)

    write(c(sum(unlist(lapply(pops, n))), m), file = filename, n = 2, sep = "\t")
    write(markers, file=filename, n = 1, append = TRUE)

    for (i in 1:length(pops)) {
      write(c("POP", c(pops[[i]]$Name)), file = filename, n = 2, sep = "\t", append = TRUE)
      write.table(twoColScoresAlt(pops[[i]], markers), file = filename, append = TRUE, quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = FALSE)
    }
  } else {

  # function for writing spam files
  writeSPAM <- function(popName, markers, counts, filename, index) {

    if (!is.null(findMarkers)) {
      for (i in 1:length(findMarkers)){
        markers <-  gsub(findMarkers[i], replaceWith[i], markers)
      }
    }
    if (!is.null(removePrefix)) markers <- gsub(removePrefix, "", markers)
    write(paste("# ", index, " ", popName, sep = ""), file=filename, append=index!=1)
    rownames(counts) <- format(substr(markers, 1, 10), width = 20, justify = "left")

    counts <- apply(counts, 2, format, width = 4, justify = "right")
    write.table(counts[,2:3], file=filename, append=TRUE, sep =" ", col.names = FALSE, quote = FALSE)
  }

# function for writing bayes files
  writeBAYES <- function(popName, markers, counts, filename, index) {

    outTable <- cbind(index, 1:length(markers), counts)
    outTable <- apply(outTable, 2, format, width = 4, justify = "right")
    write.table(outTable, file=filename, append=index!=1, sep =" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }



## iterate through populations for bayes and spam
  for (i in 1:length(pops)) {

  # setup counts table
    counts <- array("0", dim = c(m, 3))

  # iterate through markers
    for (j in 1:m) {

    # fill counts table
      countA1 <- sum(pops[[i]]$Scores[,markers[j],] == markerAlleles[markers[j],1])
      countA2 <- sum(pops[[i]]$Scores[,markers[j],] == markerAlleles[markers[j],2])

      counts[j,1] <- countA1 + countA2
      counts[j,2] <- countA1
      counts[j,3] <- countA2

    }# marker loop

  # write to file for population
  # format markers for SPAM (10 characters max)



    if (fileType == "SPAM")  writeSPAM(pops[[i]]$Name, markers, counts, filename, i)
    if (fileType == "BAYES") writeBAYES(pops[[i]]$Name, markers, counts, filename, i)

  } # population loop
}  # else end

print(paste("A", fileType, "baseline file has been written to", paste(getwd(), filename, sep ="/")))


}


