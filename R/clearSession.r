#' @title Clear \pkg{idfgen}
#'
#' @description This function clears the R console of all things \pkg{idfgen} related
#'
#' @author Eric Grau and Mike Ackerman
#'
#' @export
#' @return NULL

clearSession <- function() {
#######################################3
#
# This function clears the R console of all things IDFGEN related

  rm(list=ls(pos = 1, all=TRUE), pos=1)

  if("IDFGEN_Data" %in% search()) detach(IDFGEN_Data)
  if("Pops" %in% search()) detach(Pops)
  if("IDFGEN" %in% search()) detach(IDFGEN)



}
