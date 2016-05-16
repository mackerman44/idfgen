#' @title Resolve "AB"-"BA" discrepancies
#'
#' @description \code{readInData()} the \code{readInData()} function checks for any "AB"-"BA" discrepancies during genotype import.
#' The \code{basepairSwitches()} function can be used to convert any "BA" genotypes to "AB" genotypes
#'
#' @param x use default \code{x = Populations}
#' @param filename what would you like to name the export summary
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL


basepairSwitches <- function(x, filename=NULL) {


  allscores <- scores(x, popCol = TRUE)

# list of observed genotype counts (remove pop column is the -1)
  genotypes <- apply(allscores[,-1], 2, table, exclude = "00")

# coersion of data types here: genotypes could be array or list.  arrays will conatin no switches
  counts <- lapply(genotypes, length)

  if (4 %in% counts) {
    switchedLoci <- names(counts[counts==4])
    gen <- genotypes[switchedLoci]

  # find populations with switch - strealined function, not the prettiest

    findSwitch <- function(g, a) {  # x will be pulled from .GlobalEnv (not good programming)
      a <- names(a)
      l <- list(unique(x[x[,g]==a[2],1]), unique(x[x[,g]==a[3],1]))
      names(l) <- a[2:3]
      l
    }
    x<-allscores  # for simplicity

    switchedPops <- mapply(findSwitch, as.list(switchedLoci), gen, SIMPLIFY=FALSE)
    names(switchedPops) <- names(gen)

    if(!is.null(filename)) {

    x <- switchedPops
    rows <- max(unlist(lapply(x, lapply, length)))
    fill <- function(x, n) c(x, rep(NA, n - length(x)))
    a <- lapply(x, lapply, fill, rows)
    d <- function(x) do.call(cbind, x)
    out <- d(lapply(a, d))
    out[is.na(out)] <- ""
    write(paste(names(x), "\t", sep = ""),file=filename,sep = "\t", ncolumns=length(x)*2-1)
    suppressWarnings(write.table(out, file=filename,col.names=TRUE,row.names=FALSE, quote=FALSE,sep="\t", , append=TRUE))
    }

    return(switchedPops)

  }
  else {
    return("No switches found.")
  }


}

