#' @title Export GSI mixture file
#'
#' @description \code{dumpMixture()} exports a gsi mixture file in BAYES, SPAM, or gsi_sim format
#'
#' @param x the 'population' to be included in the gsi mixture file
#' @param markers the markers to be included in gsi mixture file. NOTE: should be same markers set that was used for the \code{dumpBaselin()}
#' function
#' @param fileType the format of the exported gsi mixture file. Can be set to "SPAM", "BAYES", or "gsi_sim"
#' @param filename what would you like to name the exported gsi mixture file
#' @param individuals keep as default \code{NULL}
#' @param findMarkers i don't remember what this does, just leave as the default \code{NULL}
#' @param replaceWith something related to \code{findMarkers}. Also leave as the default \code{NULL}
#' @param removePrefix this can be used to trim the marker names for "SPAM"
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpMixture()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpMixture <- function(x, markers, fileType, filename, individuals = NULL, findMarkers = NULL, replaceWith = NULL, removePrefix = NULL, replaceBP = TRUE, basePairs = c("A", "C", "G", "T", "-"), replacements = c("1", "2", "3", "4", "5")){
###############################################################################
#
#   Outputs mixture in gsi_sim, BAYES or SPAM format
#
#    updated 12/29/11 to include gsi_sim export ability
#    updated 2/21/12 to for version 2.1 by EG


# set individuals to export
  if(is.null(individuals)) {
    inds <- inds(x)
  } else {
    inds <- individuals
  }
  n <- length(inds)
  m <- length(markers)

  if (fileType == "gsi_sim") {


    write(c(n, m), f = filename, n = 2, a = FALSE, sep = "\t")
    write(markers, f = filename, n = 1, a = TRUE, sep = "\t")
    write(c("POP", x$Name),f = filename, n = 2, a = TRUE, sep = "\t")

    if(replaceBP) x <- replaceBPs(x, basePairs, replacements)

    output <- cbind(array(paste(x$Scores[inds,markers,1], x$Scores[inds,markers,2], sep = " "), dim = c(n, m)))
    rownames(output) <- inds

    write.table(output, file = filename, append = TRUE, quote = FALSE, sep = "    ",
              row.names = TRUE,
              col.names = FALSE)



  } else {

  #  setup markerCounts
    markerCount <- array("00", dim = c(n, m))
    rownames(markerCount) <- inds
    colnames(markerCount) <- markers

  # iterate through markers
   for(k in 1:length(markers)){

  # iterate through individuals
    for (j in 1:length(inds)){

    # get and sort individual's alleles for a given marker
      alleles <- sort(c(x$Scores[inds[j], markers[k] , 1], x$Scores[inds[j], markers[k] , 2]))

    # fill markerCounts table with genotype
      markerCount[j,k] <- paste(alleles[1], alleles[2], sep = "")

    # get marker's possible genotypes
      homo1 <- paste(markerAlleles[markers[k],1], markerAlleles[markers[k],1], sep = "")
      het <- paste(markerAlleles[markers[k],1], markerAlleles[markers[k],2], sep = "")
      homo2 <- paste(markerAlleles[markers[k],2], markerAlleles[markers[k],2], sep = "")

    # replace genotypes with counts
      markerCount[,markers[k]][markerCount[,markers[k]]==homo1] <- "20"
      markerCount[,markers[k]][markerCount[,markers[k]]==het] <- "11"
      markerCount[,markers[k]][markerCount[,markers[k]]==homo2] <- "02"
    } # indiviual loop
  } # marker loop

# write files
  if(fileType == "BAYES") write.table(markerCount, row.names = FALSE, col.names = FALSE, quote = FALSE, file = filename, sep = " ")

  if(fileType == "SPAM"){
    write("* characters", file=filename)
    if (!is.null(findMarkers)) {
      for (i in 1:length(findMarkers)){
        markers <-  gsub(findMarkers[i], replaceWith[i], markers)
      }
    }

    if (!is.null(removePrefix)) markers <- gsub(removePrefix, "", markers)
    write(paste(1:m, substr(markers, 1, 10)), file = filename, append = TRUE)
    write("* end", file = filename, append = TRUE)
    write("\\" , file = filename, append = TRUE)
    write.table(markerCount, row.names = FALSE, col.names = FALSE, quote = FALSE, file = filename, sep = " ", append = TRUE)
  }

} # end BAYES or SPAM

  print(paste("A", fileType, "mixture file has been written to", filename))
#end function

}


