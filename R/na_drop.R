#' Remove NA values
#'
#' Removes missing values from a dataset.
#'
#' @param xin A numeric vector.
#' @returns A vector without missing values
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#' @export
#' @keywords internal
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' na.drop(x)
#'
na.drop<-function(xin){
  noNA<-as.numeric(length(which(is.na(xin)==T)))
  if(noNA>0){
    x<-as.numeric(na.omit(xin))
    message(paste(noNA," missing value(s) removed"))
  }else{
    x<-xin
  }
  return(x)
}
