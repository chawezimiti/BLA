#' Summary plots
#'
#' A function to produce summary plots of a set of data.
#'
#' @param x A vector of numeric values.
#' @param varname The name of the variable (optional), character so in quotes e.g.
#'   "Clay content". If not used then the variable is called x on plots.
#' @returns A histogram with a boxplot over it and QQ plot of the variable x.
#' @author Richard Murray Lark <murray.lark@@nottingham.ac.uk>
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' summaplot(x)
#' @export
#'
summaplot<-function(x,varname){

  # plot a histogram with boxplot and QQ plot of data in x indicating
  # any probable outliers by Tukey's criterion

  x<-na.drop(x)
  if(missing(varname)){varname<-"x"}

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
  par(mfrow=c(1,2))
  ymax<-max((hist(x,plot=F))$counts)
  hist(x,main="",col="AliceBlue", xlab=varname,ylim=c(0,(ymax*1.25)))

  boxmin<-ymax*1.1
  boxmax<-ymax*1.2
  boxmid<-ymax*1.15

  lines(c(Q1,Q3),c(boxmin,boxmin))
  lines(c(Q1,Q3),c(boxmax,boxmax))
  lines(c(Q1,Q1),c(boxmin,boxmax))
  lines(c(Q3,Q3),c(boxmin,boxmax))
  lines(c(Q1,lw),c(boxmid,boxmid))
  lines(c(Q3,uw),c(boxmid,boxmid))
  lines(c(Q2,Q2),c(boxmin,boxmax),lwd=2)

  lines(c(Fu,Fu),c(10,boxmid),lty=5,col="red")
  lines(c(Fl,Fl),c(10,boxmid),lty=5,col="red")

  qqn<-qqnorm(x,main="",pch=16)
  qqline(x)
  points(qqn$x[ol],qqn$y[ol],pch=16,col="red")

  par(mfrow=c(1,1))
}
