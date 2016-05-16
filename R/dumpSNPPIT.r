#' @title Export snppit file
#'
#' @description \code{dumpSNPPIT()} exports a file in SNPPIT format
#'
#' @param baselinePops the list of populations to be used as the parental baseline
#' @param mixturePops the populations to be included as potential offspring
#' @param markers the list of markers to be used in the snppit file
#' @param filename what would you like to name the exported snppit file?
#' @param errorFile if you have a .txt file containing the genotype error rates for each of the markers to be used in the snppit file, what is the
#' name of the file? The .txt should contain 2 columns: 1) the name of each markers and 2) the error rate to be used
#' @param errorDefaultValue if you don't have error rates for each of your markers, what would you like snppit to use as the default error
#' rate
#' @param POPCOLUMN_SEX If you would like to include sex information for the parents, which column in your input file contains the sex data?
#' @param POPCOLUMN_REPRO_YEARS If you would like to include the year that the parents reproduced, which column in your input file contains the
#' reporductive years?
#' @param POPCOLUMN_SPAWN_GROUP If you would like to include spawn groups for the parents, which column in your input file contains the spawn groups?
#' @param OFFSPRINGCOLUMN_BORN_YEAR If you would like to include the year born for offspring, which column contains that information?
#' @param OFFSRPINGCOLUMN_SAMPLE_YEAR If you would like to include sampling year for offspring, which column in your input file contains that
#' information?
#' @param OFFSPRINGCOLUMN_AGE_AT_SAMPLING If you would like to include the age at sampling for offspring, which column in your input file
#' contains that information?
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpSNPPIT <- function(baselinePops,
                      mixturePops,
                      markers,
                      filename,
                      errorFile = NULL,
                      errorDefaultValue = NULL,
                      POPCOLUMN_SEX = NULL,
                      POPCOLUMN_REPRO_YEARS = NULL,
                      POPCOLUMN_SPAWN_GROUP = NULL,
                      OFFSPRINGCOLUMN_BORN_YEAR = NULL,
                      OFFSRPINGCOLUMN_SAMPLE_YEAR = NULL,
                      OFFSPRINGCOLUMN_AGE_AT_SAMPLING = NULL) {


# build error table

# CASE 1  Error Table with defaults if needed
  if (!is.null(errorFile)) {
  # get error table and format
    errors <- read.table(errorFile, sep = "\t", colClasses = "character", row.names = 1)
    rownames(errors) <- gsub("[.]", "",  rownames(errors))
    rownames(errors) <- gsub("[-]", "",  rownames(errors))

  # compare table with markerList
    inTableErrors <- intersect(rownames(errors), markers)
    notInTableErrors <- setdiff(markers, rownames(errors))

  # use default values if needed
    if (length(notInTableErrors > 0)) {
    # print markers not in error table
      print(notInTableErrors)
      print("Warning! The above markers are not in the error table.")

    # check for default value
      if (!is.null(errorDefaultValue)) {
      # use default value
        print(paste("The default value of", errorDefaultValue, "will be used."))
        tempErrorTable <- c(errors[inTableErrors,1], rep(errorDefaultValue, length(notInTableErrors)))
        errorTable <- array(tempErrorTable,  dim = c(length(tempErrorTable), 1), dimnames = list(c(inTableErrors, notInTableErrors), "Error"))

      } else {
      # quit if no value and send warning
        return (print("WARNING!!! FUNCTION EXIT!!! Please provide a errorDefaultValue in the function call."))
      }
  # no default values needed
    } else {
     errorTable <- array(errors[inTableErrors,1],  dim = c(length(inTableErrors), 1), dimnames = list(inTableErrors, "Error"))

    }
# CASE 2 No error file
  } else if (!is.null(errorDefaultValue)) {
    print(paste("No error table included.  The default value of", errorDefaultValue, "will be used"))
    errorTable <- array(as.character(errorDefaultValue), dim = c(length(markers), 1), dimnames = list(markers, "Error"))
# CASE 3 NO BOTH - QUIT with warning
  } else {

    return (print("You must provide the argument errorFile or defaultErrorValue or both."))
  }


# write header
  write(paste("NUMLOCI", length(markers)), filename)
  write("MISSING_ALLELE 0", filename, append = TRUE)
  if(!is.null(POPCOLUMN_SEX)) write("POPCOLUMN_SEX", filename, append = TRUE)
  if(!is.null(POPCOLUMN_REPRO_YEARS)) write("POPCOLUMN_REPRO_YEARS", filename, append = TRUE)
  if(!is.null(POPCOLUMN_SPAWN_GROUP)) write("POPCOLUMN_SPAWN_GROUP", filename, append = TRUE)
  if(!is.null(OFFSPRINGCOLUMN_BORN_YEAR)) write("OFFSPRINGCOLUMN_BORN_YEAR", filename, append = TRUE)
  if(!is.null(OFFSRPINGCOLUMN_SAMPLE_YEAR)) write("OFFSRPINGCOLUMN_SAMPLE_YEAR", filename, append = TRUE)
  if(!is.null(OFFSPRINGCOLUMN_AGE_AT_SAMPLING)) write("OFFSPRINGCOLUMN_AGE_AT_SAMPLING", filename, append = TRUE)




# write error table
  errorTable <- errorTable[markers,]
  write.table(errorTable, file = filename, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE)

# write populations function
  writeInds <- function(popNames, pops = TRUE) {


    for (i in 1:length(popNames)) {
    # format data


      popData <- get(popNames[i])
      popData$Scores[popData$Scores=="A"] <- 1
      popData$Scores[popData$Scores=="C"] <- 2
      popData$Scores[popData$Scores=="G"] <- 3
      popData$Scores[popData$Scores=="T"] <- 4
      popData$Scores[popData$Scores=="-"] <- 5

      if (pops) {


        output <- cbind(rownames(popData$Scores),
                        if(!is.null(POPCOLUMN_SEX)) popData$IndividualData[,POPCOLUMN_SEX],
                        if(!is.null(POPCOLUMN_REPRO_YEARS)) popData$IndividualData[,POPCOLUMN_REPRO_YEARS],
                        if(!is.null(POPCOLUMN_SPAWN_GROUP)) popData$IndividualData[,POPCOLUMN_SPAWN_GROUP],
                        array(paste(popData[[1]][,markers,1], popData[[1]][,markers,2], sep = "\t"), dim = c(length(popData[[1]][,1,1]), length(markers))))
      } else {

          inds <- rownames(popData$Scores)


        if (length(inds) > 0 ) {
        output <- cbind(rownames(popData$Scores[inds,,]),
                        if(!is.null(OFFSPRINGCOLUMN_BORN_YEAR)) popData$IndividualData[inds,OFFSPRINGCOLUMN_BORN_YEAR],
                        if(!is.null(OFFSRPINGCOLUMN_SAMPLE_YEAR)) popData$IndividualData[inds,OFFSRPINGCOLUMN_SAMPLE_YEAR],
                        if(!is.null(OFFSPRINGCOLUMN_AGE_AT_SAMPLING)) popData$IndividualData[inds,OFFSPRINGCOLUMN_AGE_AT_SAMPLING],
                        array(paste(popData[[1]][inds,markers,1], popData[[1]][inds,markers,2], sep = "\t"), dim = c(length(popData[[1]][inds,1,1]), length(markers))))
        }
      }

      if (pops){
        popLine <- paste("POP", names(popData), sep = " ")
        write(popLine, filename, append = TRUE)
        write.table(output, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      } else if (length(inds) > 0) {
        popLine <- paste(paste("OFFSPRING", names(popData)), "?", sep = "\t")
        write(popLine, filename, append = TRUE)
        write.table(output, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
      }


    }
  }

  writeInds(baselinePops)
  writeInds(mixturePops, FALSE)






#errorTable
}
