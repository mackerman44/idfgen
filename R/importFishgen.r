#' @title Import Data from www.fishgen.net
#'
#' @description \code{readInData()} imports data exported from www.fishgen.net
#'
#' @param inputFile the .txt file from www.fishgen.net containing the biological and genotype data to be imported into \pkg{idfgen}.
#' @param genotypeStartColumn the column # in \code{inputFile} where the genotypes begin. This will vary depending on the number of additional fields added to the
#' Progeny export.
#' @param popColumn which column in \code{inputFile} would you like to use to parse the data. The default is 1, which is where Progeny exports Pedigree
#' names
#' @param sortBPs the \code{readInData()} function will check for "AB"-"BA" discrepancies. If \code{sortBPs = TRUE}, then any "BA" genotypes will
#' automatically be converted to "AB"
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

importFishgen <- function(inputFile, genotypeStartColumn, popColumn = 1, sortBPs = FALSE)
{
#
#  This function takes a progeny export and creates objects needed for
#  further analysis
#
#  written 01/12/2012 by EG
#
#

# read file
  rawData <- read.table(inputFile, quote="", header = TRUE, sep="\t", row.names = 3, as.is = TRUE, colClasses = "character", comment.char = "")

# adjust column positions since ind names column is removed
  genotypeStartColumn <- genotypeStartColumn - 1

# adjust column position of popColumn accordingly
  if ((length(popColumn) != 1) | !(1 %in% popColumn)) popColumn <- popColumn - 1

# fill zeroes if blanks
  rawData[,genotypeStartColumn:length(colnames(rawData))][rawData[,genotypeStartColumn:length(colnames(rawData))]==""] <- 0


# marker list
  markers <- colnames(rawData)[(genotypeStartColumn):length(colnames(rawData))]
  allele1Index <- 1:(length(markers)/2) * 2 - 1
  allele2Index <- 1:(length(markers)/2) * 2
  #markers <- gsub(".A1", "", markers)      # remove suffix for progeny export
  markers <- substr(markers,1,nchar(markers)-3)
  markers <- gsub("[.]", "", markers)      # remove decimals
  markers <- markers[allele1Index]

# setup markerAlleles table
  ma <- rawData[,genotypeStartColumn + allele1Index - 1]
  mb <- rawData[,genotypeStartColumn + allele2Index - 1]
  colnames(ma) <- markers
  colnames(mb) <- markers
  rbind(ma, mb) -> m

  maxAlleleCount <- max(unlist(lapply(apply(m, 2, unique), length)))
  if (maxAlleleCount == 1) maxAlleleCount <- 2
  markerAlleles <- array(dim = c(length(markers), maxAlleleCount))


# iterate through markers to build markerAllele table
  for (m in 1:length(markers)) {

  # get all alleles for a marker
    alleles <- unique(c(rawData[,genotypeStartColumn + (2 * m) - 2], rawData[,genotypeStartColumn + (2*m) - 1]))
  # ensure consistency by sorting
    alleles <- sort(alleles)

  # remove 0 if present
    alleles <- alleles[alleles!=0]

  # assign alleles to table
    if (length(alleles) > 0) {
      markerAlleles[m,1:length(alleles)] <- alleles

      if (length(alleles) < 2) {
        markerAlleles[m, 2] <- "monomorphic"
      }
    }
  # update table
  } # marker loop

# finish markerAlleles setup
  rownames(markerAlleles) <- markers
  colnames(markerAlleles) <- paste("Allele", 1:maxAlleleCount, sep = "")
  noDataMarkers <- rownames(markerAlleles)[is.na(markerAlleles[,1])==TRUE]
  if (length(noDataMarkers) > 0) {
    cat("**** The following markers have no data associated with them ****\n")
    print(noDataMarkers)
  }

  markerLength <- nchar(markerAlleles[,1])
  uSatMarkers <- names(markerLength[markerLength > 1])
  snpMarkers <- markers[!(markers %in% uSatMarkers)]
  uSatMarkers <- uSatMarkers[!(uSatMarkers %in% noDataMarkers)]

# population names

  if (length(popColumn) == 1) {

    # format for R variable name
    popNamesForAllInds <- make.names(rawData[,popColumn])
    popColumnName <- make.names(colnames(rawData)[popColumn])

    # put column in first position
    rawData <- cbind(as.character(popNamesForAllInds), rawData[,-popColumn])
    rawData[,1] <- as.character(rawData[,1])
    colnames(rawData)[1] <- popColumnName
  } else {


    # combine columns and format for R
    popNamesForAllInds <- make.names(paste(rawData[,popColumn[1]], rawData[,popColumn[2]]))
    popColumnName <- make.names(paste(colnames(rawData)[popColumn[1]], colnames(rawData)[popColumn[2]]))

    # add column to first position
    rawData <- cbind(popNamesForAllInds, rawData)
    rawData[,1] <- as.character(rawData[,1])
    colnames(rawData)[1] <- popColumnName

    # increase column count to account for additional column
    genotypeStartColumn <- genotypeStartColumn + 1

  }

  popNames <- sort(unique(rawData[,1]))



# metaData
  metaDataFields <- colnames(rawData)[1:(genotypeStartColumn - 1)]


# For each population create object

  allPops <- list()
  for (i in 1:length(popNames)){

  # partition for population i
    popData <- rawData[rawData[,1]==popNames[i],]

  # get genotype data
    allele1 <- popData[allele1Index + (genotypeStartColumn - 1)]
    allele2 <- popData[allele2Index + (genotypeStartColumn - 1)]
    colnames(allele1) <- markers
    colnames(allele2) <- markers

    if(sortBPs) {
      d <- allele1
      d[allele1 > allele2] <- allele2[allele1 > allele2]
      allele2[allele1 > allele2] <- allele1[allele1 > allele2]
      allele1 <- d
    }

  # create [ind,marker,allele] 3 dimensional array  (the Score)
    scores <- abind(allele1, allele2, along = 3)  # uses package abind here
    dimnames(scores) <- list(Inds = rownames(allele1), Markers = markers, Allele = c("Allele1", "Allele2"))

  # metadata
    individualData <- popData[metaDataFields]

  # info
   #info <- list(popName = popNames[i], timeStamp = Sys.Date())


  # output population.data
    tempPop <- list()

    tempPop$Scores <- scores
    tempPop$IndividualData <- individualData
    tempPop$Name <- popNames[i]
    tempPop$Individuals <- rownames(scores)
    tempPop$Markers <- colnames(scores)

    class(tempPop) <- "Population"
    allPops[[i]] <- tempPop
    #assign(popNames[i], tempPop, pos = 1)



  } # population loop

  populations <- popNames
  class(populations) <- "PopList"
  names(allPops) <- popNames
  attach(allPops, name = "IDFGEN_Data")

# accesssory outputs

 # IDFGEN_Data <- new.env()

  assign("markerAlleles", markerAlleles, pos = "IDFGEN_Data")
  #attach(IDFGEN_Data)
 # rm(IDFGEN_Data)

  #assign("markerListAll", markers, pos = 1)
  #assign("popNamesAll", popNames, pos = 1)

  assign("MetaDataFields", metaDataFields, pos = 1)
  assign("MarkersSNP", snpMarkers, pos = 1)
  assign("MarkersUSAT", uSatMarkers, pos = 1)

  assign("Markers", markers, pos = 1)
  assign("Populations", populations, pos = 1)

# print to console
  cat(paste("****", length(popNames), "new populations created!", "****\n"))
  print(popNames, ncol = length(popNames))
  cat("\n**** Objects created ****\n")
  cat("Populations    : a 'PopList' of the names of all populations\n")
  cat("Markers        : a vector of the names of all markers\n")
  cat("MarkersSNP     : a vector of the names of all SNPs\n")
  cat("MarkersUSAT    : a vector of the names of all uSATs\n")
  cat("MetaDataFields : a vector of the names of all meta data fields\n\n")


# run basepair switches
if (sortBPs == FALSE) {

  cat("Checking for switched base pairs...")
  bpswitch <- basepairSwitches(Populations)
  if(class(bpswitch)=="list") {
    cat("\n** Warning:", length(bpswitch), "markers were found with 'AB' - 'BA' allele switches\n\n")
    cat("** Use basepairSwitches(Populations, filename ='myfile.txt') to export results.\n\n")
    cat("** Use readInData(..., sortBPs=TRUE) to ensure consistent allele order.\n")
    cat("** Recommended: run clearSession() to remove all objects, populations, packages, etc and start the project analysis over if consistent allele order is desired..\n\n")
  }
  else cat("no errors found.\n")

}

}
