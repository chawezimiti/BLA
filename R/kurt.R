#' kurtosis
#'
#' A function to calculate an estimate of the coefficient of kurtosis
#' from a set of data.
#'
#' @param x A vector of numeric values.
#'
#' @returns The reduced coefficient of kurtosis.
#'
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#'
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' kurt(x)
#'
#' @keywords internal
#' @export
#'
kurt<-function(x){

  x<-na.drop(x)
  n<-length(x)
  xd<-x-mean(x)
  mu4<-sum(xd^4)/(n-1)
  mu<-sqrt(sum(xd^2)/(n-1))
  sk<-(mu4/(mu^4))-3
  return(sk)
}
