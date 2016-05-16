#' @title Write a simple table
#'
#' @description \code{dumpTable()} just exports any simple table, perhaps containing the biological data for your individuals
#'
#' @param x the object to be written, preferably a matrix or data frame. If not, it is attempted to coerce \code{x} to a data frame.
#' @param filename either a character string naming a file or a \link{connection} open for writing. "" indicates output to the console.
#' @param row.names either a logical value indicating whether the row names of \code{x} are to be written along with \code{x}, or a character
#' vector of row names to be written
#' @param sep the field separator string. Values within each row of x are separated by this string.
#'
#' @seealso \code{\link[utils]{write.table}}
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

dumpTable <- function(x, filename, row.names = FALSE, sep = "\t") {

  write.table(x, file = filename, append = FALSE, quote = FALSE, sep = sep,
            row.names = row.names, col.names = if(row.names) NA else TRUE)


}
