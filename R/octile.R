#' Octile Skewness
#'
#' A function to calculate an estimate of the octile skewness from a set of data.
#'
#' @param x A vector of numeric values.
#' @returns The octile skewness.
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#' @keywords internal
#' @export
#'
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' ocskew(x)
#'
ocskew<-function(x){

  x<-na.drop(x)
  Ocs<-quantile(x,c(1/8,0.5,7/8))
  os<-((Ocs[3]-Ocs[2])-(Ocs[2]-Ocs[1]))/(Ocs[3]-Ocs[1])
  return(os)
}
