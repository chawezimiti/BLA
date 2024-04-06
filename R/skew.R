#' Skewness
#'
#' Function to calculate an estimate of the coefficient of skewness from
#' a set of data.
#'
#' @param x A vector of numeric values.
#' @returns The coefficient of skewness.
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#' @keywords internal
#' @export
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' skew(x)
#'
skew<-function(x){

  x<-na.drop(x)
  n<-length(x)
  xd<-x-mean(x)
  mu3<-sum(xd^3)/(n-1)
  mu<-sqrt(sum(xd^2)/(n-1))
  sk<-mu3/(mu^3)
  return(sk)
}

