#####
# generics.r
#
#

# CUSTOM FUNCTION INITIALIZERS


# Accessors
  n <- function(x, ...) UseMethod("n")
  m <- function(x, ...) UseMethod("m")
  inds <- function(x, ...) UseMethod("inds")
  markers <- function(x, ...) UseMethod("markers")
  metaData <- function(x, ...) UseMethod("metaData")
  oneColScores <- function(x, ...) UseMethod("oneColScores")
  twoColScores <- function(x, ...) UseMethod("twoColScores")
  twoColScoresAlt <- function(x, ...) UseMethod("twoColScoresAlt")
  scores <- function(x, ...) UseMethod("scores")


# Data Manipulation
  findDuplicateInds <- function(x, ...) UseMethod("findDuplicateInds")
  findNoCalls <- function(x, ...) UseMethod("findNoCalls")
  findAlleles <- function(x, ...) UseMethod("findAlleles")
  removeIndividuals <- function(x, ...) UseMethod("removeIndividuals")
  poolPops <- function(x, ...) UseMethod("poolPops")
  subPopulation <- function(x, ...) UseMethod("subPopulation")
  replaceBPs <- function(x, ...) UseMethod("replaceBPs")

# Data Output

  



# Default & Override Functions
  removeIndividuals.default <- function(x, ...) {
    if((class(x)=="matrix") | (class(x)=="array") | (class(x)=="data.frame"))  {
      pops <- unique(x[,1])
      total <- 0
      for (i in 1:length(pops)) {
        nRemoved <- removeIndividuals(get(pops[i]), x[,2][x[,1]==pops[i]])
        total <- total + nRemoved
      }
      cat("\n***", total, "total individuals have been removed *** \n")
    } else {
      cat("***", deparse(match.call()[-1]), "is an invalid data type -- see documentation\n")
    }   
  }

findHybrids <- function(...) {
  findAlleles(...)
}



#### Population Genetic Analyses


npsScores <- function(x,...) UseMethod("npsScores")
#maf <- function(x, ...) UseMethod("maf")












