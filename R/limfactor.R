#' Determination of the most limiting factor to biological response
#'
#' This function determines the most limiting factor based on von Liebig law of
#' the minimum given results of the predicted boundary line values for the different
#' factors of interest. Boundary lines for various factors are fitted and the factor
#' that predicts the minimum response for a particular point is considered as the
#' most limiting factor (Casanova et al. 1995).
#'
#' @param ... vectors with predicted values from the boundary line models for
#'   each factor being evaluated.
#' @returns A dataframe consisting of the most limiting factor and the minimum
#'   predicted response
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#' @export
#' @examples
#'
#' N<-rnorm(10,50,5)#assuming these are predicted responses using the fitted BL for N,P,K
#' K<-rnorm(10,50,4)
#' P<-rnorm(10,50,6)
#'
#' limfactor(N,K,P)
#'
limfactor<-function(...){

  dat<-data.frame(...)
  name<-names(dat)
  dat2<- as.data.frame(t(dat)) # transpose the dataframe. changes rows to columns and vice-versal
  dat3<-suppressWarnings(unlist(lapply(dat2,min,na.rm=T))) # minimum values for each column
  dat3[which(dat3==Inf)]<-NA # Inf is returned if the column has only NA's.
  dat4<-lapply(dat2,which.min)

  dat5<-lapply(dat4, function(a,b){
    out<-b[a]
    return(out)
  },b=name)

  dat6<-lapply(dat5,function(b){ifelse(length(b)==0,b<-NA,b<-b)})

  data<-data.frame(Rs=dat3,Factor=unlist(dat6))
  data2<-data.frame(Rs=data$Rs, Lim_factor=data$Factor)

  ## Determining unidentified yield gap

  test<-apply(dat2,1,max, na.rm=T)

  unidentified<-unlist(lapply(dat2,function(x,test){
    all(x==test)
  },test))

  for(i in 1:dim(data2)[1]){
    ifelse(unidentified[i]==TRUE,data2$Lim_factor[i]<-"unidentified",data2$Lim_factor[i]<-data2$Lim_factor[i])
  }

  ##

  Largest<-max(dat, na.rm=T)

  return(list(data2,Largest))

}
