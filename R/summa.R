#' Summary statistics
#'
#' A function to calculate summary statistics of a set of data.
#'
#' @param x A vector of numeric values.
#' @param sigf The number of significant figures to report (optional).
#' @returns A matrix containing the mean value, median value,
#'   first and third quartiles, sample variance, sample standard deviation,
#'   coefficient of skewness, octile skewness, coefficient of kurtosis and
#'   the number of probable outliers in a data set.
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' summa(x,2)
#' @export
#'
summa<-function(x,sigf){

  # compute summary statistics of values in x

  if(missing(sigf)){rosig<-F}else{rosig<-T}


  x<-na.drop(x)
  Q1<-quantile(x,prob=0.25)
  Q3<-quantile(x,prob=0.75)
  Q2<-quantile(x,prob=0.5)
  hspread<-Q3-Q1
  Fu<-Q3+3*hspread
  Fl<-Q1-3*hspread

  # ols,oll: values below and above outer fences
  # posols,posoll: values below and above inner fences
  # (so ols and posols overlap, as do oll and posoll
  #

  ols<-which(x<Fl)
  oll<-which(x>Fu)
  posols<-which(x<(Q1-1.5*hspread))
  if(length(posols)==0){
    lw<-min(x)}else{
      lw<-min(x[-posols])}
  posoll<-which(x>(Q3+1.5*hspread))
  if(length(posoll)==0){
    uw<-max(x)}else{
      uw<-max(x[-posoll])}

  ol<-c(ols,oll) # combined outlier set

  nol<-length(ol)

  outp<-matrix(c(mean(x),Q2,Q1,Q3,var(x),sqrt(var(x)),skew(x),ocskew(x),kurt(x),nol),1,10)

  if(rosig=="TRUE"){outp<-signif(outp,sigf)}

  colnames(outp)<-c("Mean","Median",
                    "Quartile.1", "Quartile.3","Variance","SD","Skewness",
                    "Octile skewness","Kurtosis",
                    "No. outliers")


  return(outp)

}
