#' @title Export COLONY file
#'
#' @description \code{dumpColony()} exports a file in Colony format
#'
#' @param x the list of populations to be included in Colony file
#' @param markers the markers to be included in Colony file
#' @param filename what would you like to name the exported COLONY file?
#' @param errorFile what would you like to name the marker error rate file for COLONY?
#' @param errorDefaultValues what would you like your defaule genotype error rate to be for COLONY?
#' @param replaceBP would you like to replace the allele names with something else. By default, \code{dumpColony()} replaces
#' c("A","C","G","T","-") with c("1","2","3","4","5"). This can be set to "FALSE" to keep data in the original imported format. Alternatively,
#' the \code{basePairs} and \code{replacements} arguments can be modified to accommodate any desired import or export format.
#' @param basePairs If \code{replaceBP = TRUE}, the nomenclature used for the imported genotype data
#' @param replacements If \code{replaceBP = TRUE}, what would you like to replace the alleles with?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpColony <- function(x,
                       markers = NULL,
                       filename,
                       errorFile = NULL,
                       errorDefaultValues = c(0,0),
                       replaceBP = TRUE,
                       basePairs = c("A", "C", "G", "T", "-"),
                       replacements = c("1", "2", "3", "4", "5"))
{


  errorTable <- array(NA, dim= c(length(markers), 2), dimnames = list(markers, c("e1", "e2")))

# FILE VALUES + DEFAULTS if needed
  if (!(is.null(errorFile))) {
    errors <- read.table(errorFile, sep = "\t", colClasses = "character", row.names = 1)
    rownames(errors) <- gsub("[.]", "",  rownames(errors))
    rownames(errors) <- gsub("[-]", "",  rownames(errors))
    inTableMarkers <- intersect(rownames(errors), markers)
    notInTableMarkers <- setdiff(markers, rownames(errors))
    if (length(inTableMarkers) > 0) {
      errorTable[inTableMarkers,1] <- errors[inTableMarkers,1]
      errorTable[inTableMarkers,2] <- errors[inTableMarkers,2]
    }
    if(length(notInTableMarkers) > 0) {
      errorTable[notInTableMarkers, ] <- rep(errorDefaultValues, each = length(notInTableMarkers))
    }
  }
# DEFAULTS VALUES ONLY
    else {
    errorTable[,] <- rep(errorDefaultValues, each = length(markers))

  }

# WRITE ERROR FILE
  errorFile <- paste(sub(".txt", "", filename), "_marker_error_rates.txt", sep = "")
  write(markers, errorFile, ncol = length(markers), sep = "\t")
  write(rep(0, length(markers)), errorFile, ncol = length(markers), append = TRUE, sep = "\t")
  write(sprintf("%1.4f", as.numeric(errorTable[,1])), errorFile, ncol = length(markers), append = TRUE, sep = "\t")
  write(sprintf("%1.4f", as.numeric(errorTable[,2])), errorFile, ncol = length(markers), append = TRUE, sep = "\t")


#########

# Write Genetic Data File

#print(x)

   if (replaceBP)  x <- replaceBPs(x, basePairs, replacements)


#  print(x$Scores)


 # scores <- cbind(inds(x), array(paste(x$Scores[,markers,1], x$Scores[,markers,2], sep = "\t"), dim = n(x), length(markers)))

#  print(scores)

 scores <- twoColScoresAlt(x,markers)

 write.table(as.array(scores), file = filename, append = FALSE, quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = FALSE)

 cat("A Colony file has been written to", filename, "\n")



}

######### dumpDat Function
dumpDat <- function (x, markers = NULL, filename, title,
                     replaceBP = TRUE,
                     basePairs = c("A", "C", "G", "T", "-"),
                     replacements = c("1", "2", "3", "4", "5"))
{
  #errorTable <- array(NA, dim = c(length(markers), 2), dimnames = list(markers,c("e1", "e2")))
  #if (!(is.null(errorFile))) {
  #  errors <- read.table(errorFile, sep = "\t", colClasses = "character",row.names = 1)
  #  rownames(errors) <- gsub("[.]", "", rownames(errors))
  #  rownames(errors) <- gsub("[-]", "", rownames(errors))
  #  inTableMarkers   <- intersect(rownames(errors), markers)
  #  notInTableMarkers<- setdiff(markers, rownames(errors))
  #  if (length(inTableMarkers) > 0) {
  #    errorTable[inTableMarkers, 1] <- errors[inTableMarkers,1]
  #    errorTable[inTableMarkers, 2] <- errors[inTableMarkers,2]
  #  }
  #  if (length(notInTableMarkers) > 0) {
  #    errorTable[notInTableMarkers, ] <- rep(errorDefaultValues,each = length(notInTableMarkers))
  #  }
  #}
  #else {
  # errorTable[, ] <- rep(errorDefaultValues, each = length(markers))
  #}
  #errorFile <- paste(sub(".txt", "", filename), "_marker_error_rates.txt",sep = "")
  #write(markers, errorFile, ncol = length(markers), sep = "\t")
  #write(rep(0, length(markers)), errorFile, ncol = length(markers),append = TRUE, sep = "\t")
  #write(sprintf("%1.4f", as.numeric(errorTable[, 1])), errorFile,ncol = length(markers), append = TRUE, sep = "\t")
  #write(sprintf("%1.4f", as.numeric(errorTable[, 2])), errorFile,ncol = length(markers), append = TRUE, sep = "\t")
  if (replaceBP)
    x <- replaceBPs(x, basePairs, replacements)
  scores <- twoColScoresAlt(x, markers)
  #filetemp <- unlist(strsplit(filename,split = ".dat"))
  write(paste("'","C:\\GenSoftware\\Colony\\",title,"\\",title,"'",sep=""), file = filename, append = FALSE, sep = "")
  write(paste("'","C:\\GenSoftware\\Colony\\",title,"\\",title,"'",sep=""), file = filename, append = TRUE, sep = "")
  write(paste(n(x),"! Number of offspring in the sample",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(length(markers),"! Number of loci",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(1234,"! Seed for random number generator",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=Not updating/updating allele frequency",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(2,"! 2/1=Dioecious/Monoecious species",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=Inbreeding absent/present",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=Diploid species/HaploDiploid species",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste("1  1","! 0/1=Polygamy/Monogamy for males & females",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1/2=no prior/sibship size prior/sibship complexity prior",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=Unknown/Known population allele frequency",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(1,"! Number of runs",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(2,"! 1/2/3/4 = Short/Medium/Long/VeryLong run",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=Monitor method by Iterate#/Time in second",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(10000,"! Monitor interval in Iterate# / in seconds",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(0,"! 0/1=DOS/Windows version",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(1,"! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write(paste(1,"! 0/1/2/3=Low/Medium/High/VeryHigh precision",sep = "\t"), file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(markers, file = filename, append = TRUE, sep = ",", ncol = length(markers))
  write(rep(0,length(markers)), file = filename, append = TRUE, sep = ",", ncol = length(markers))
  write(rep(format(0.0000,nsmall=4),length(markers)), file = filename, append = TRUE, sep = ",", ncol = length(markers))
  write(rep("0.001",length(markers)), file = filename, append = TRUE, sep = ",", ncol = length(markers))
  write("", file = filename, append = TRUE)
  write.table(as.array(scores), file = filename, append = TRUE,quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
  write("", file = filename, append = TRUE)
  write("0.0  0.0", file = filename, append = TRUE, sep = "\t")
  write("0  0", file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(0, file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write("", file = filename, append = TRUE)
  write("!____________NOTE____________!",file = filename, append = TRUE, sep = "\t")
  write("", file = filename, append = TRUE)
  write(paste("Input file",filename,"generated by dumpDat function in IDFGEN package created by MA on 2/21/14"),file = filename, append = TRUE, sep = "")
  write(paste("Generated by",Sys.getenv("USERNAME"),"on",Sys.time()),file = filename, append = TRUE, sep = "")
  cat("A .dat file has been written to", filename, "\n")
}

