#' Testing evidence of boundary existence in dataset
#'
#' This function determines the probability of having bounding effects in a scatter
#' plot of of \code{x} and \code{y} based on the clustering of points at the upper
#' edge of the scatter plot (Miti et al.2024). It tests the hypothesis of larger
#' clustering at the upper bounds of a scatter plot against a null bivariate normal
#' distribution with no bounding effect (random scatter at upper edges). It returns
#' the probability (p-value) of the observed clustering given that it a realization
#' of an unbounded bivariate normal distribution.
#'
#' @param x A numeric vector of values for the independent variable.
#' @param y A numeric vector of values for the response variable.
#' @param shells A numeric value indicating the number of boundary peels
#'   (default is 10).
#' @param simulations The number of simulations for the null bivariate normally
#'   distributed data sets used to test the hypothesis (default is 1000).
#' @param method This describes the measure of boundary points compaction. The methods
#'   include \code{"sd-enclidean"} for the euclidean distance standard deviation of the
#'   each boundary point to the center of data, \code{"Area"} for the perimeter around
#'   the boundary points and \code{"Perimeter"} for the area covering the boundary points.
#' @param alpha a relative measure of concavity of polygon if method is set to \code{"Area"}
#'   or \code{"Perimeter"}. 1 results in a relatively detailed shape, Infinity results in
#'   a convex hull. We recommend values in the range of 1 - 5.
#' @param plot If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
#'   is not part of output (default is \code{TRUE}).
#' @param ... Additional graphical parameters as with the \code{par()} function.
#'
#' @returns A dataframe containing the measures of peel compaction in the left and right
#'   sections of the data with their corresponding probability values.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#' @import MASS stats concaveman
#' @export
#'
#' @details
#' It is recommended that any outlying observations, as identified by the
#' \code{bagplot()} function of the \code{aplpack} package are removed from
#' the data. This is also implemented in the simulation step in the
#' \code{expl_boundary()} function.
#'
#' @references
#'
#' Eddy, W. F. (1982). Convex hull peeling, COMPSTAT 1982-Part I: Proceedings in
#' Computational Statistics, 42-47. Physica-Verlag, Vienna.
#'
#' Miti. c., Milne. A. E., Giller. K. E. and Lark. R. M (2024). Exploration of data
#' for analysis using boundary line methodology. Computers and Electronics in Agriculture
#' 219 (2024) 108794.
#'
#' Park, J.-S. and Oh, S.-J. (2012).A new concave hull algorithm and concaveness measure
#' for n-dimensional datasets.Journal of Information science and engineering,28(3):587–600.
#'
#' @usage
#' expl_boundary(x,y,shells=10,simulations=1000,method="sd-enclidean",alpha=1,
#'               plot=TRUE,...)
#'
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' y<-evapotranspiration$`yield(t/ha)`
#' expl_boundary(x,y,10,100) # recommendation is to set simulations to greater than 1000
#'
expl_boundary<-function(x,y,shells=10,simulations=1000,method="sd-enclidean",alpha=1,plot=TRUE,...){

  if(simulations>=1000) message("Note: This function might take longer to execute when running a large number of simulations.\n")


  ## Selection of the x_min and x_max index values----------------------------------------

  upper.peel<-function(peel){
    min.x.index<-which(peel[,2]==max(peel[,2][which(peel[,1]==min(peel[,1]))]) & peel[,1]==min(peel[,1])) # selection of min
    max.x.index<-which(peel[,2]==max(peel[,2][which(peel[,1]==max(peel[,1]))]) & peel[,1]==max(peel[,1])) # and max.x.index

    if(min.x.index<max.x.index){                                                                          # even when two values
      op<-peel[min.x.index:max.x.index,]                                                                  # occur. The highest is used
    }else{
      op<-rbind(peel[(min.x.index:nrow(peel)),],
                peel[(1:max.x.index),]
      )
    }
    return(op)
  }

  ## Removing NA'S from the data ---------------------------------------------------------

  data<- data.frame(x=x,y=y)
  test<-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE)

  if(length(test)>0){
    data1<-data[-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE),]}else{
      data1<-data
    }

  x<-data1$x
  y<-data1$y
  dat<-cbind(x,y)

  ## setting the output area -------------------------------------------------------------

  if(plot==TRUE){

    old_par <- par(no.readonly = TRUE) # Save the current graphical parameters
    on.exit(par(old_par))# Ensure the original graphical parameters are restored on exit

    plot_layout<-rbind(c(1,1,2),c(1,1,3))
    layout(plot_layout)
    plot(dat,...)}

  ## Determination of the convex hull-----------------------------------------------------

  peels<-list()
  left<-list()
  right<-list()
  n<-length(dat[,1]) #: sample size
  Sigma<-cov(dat)

  for(i in 1:shells){
    ch.index<-chull(dat) # find convex hull (index values for points)
    p1<-dat[ch.index,]   # extract peel 1
    p1_2<-p1[!duplicated(p1), ]
    p1.upper<-upper.peel(p1_2)


    p1.upperx<-data.frame(x=p1.upper[,1],y=p1.upper[,2])

    ifelse( p1.upperx[1,2]!=max(p1.upperx$y),
            peels[[i]]<-split(p1.upperx,ifelse(p1.upperx$x<p1.upperx$x[which(p1.upperx$y==max(p1.upperx$y))],0,1)),
            peels[[i]]<-split(p1.upperx,ifelse(p1.upperx$x<= p1.upperx$x[which(p1.upperx$y==max(p1.upperx$y))],0,1)))

    index.func<-function(x,y){
      del<-list()
      for(i in 1:length(y[,1])){
        del[[i]]<-which(x[,1]==y[,1][i]& x[,2]==y[,2][i])
      }
      del

      d<-unlist(del, recursive = TRUE, use.names = TRUE)

      d
    }


    dat<-dat[-index.func(dat,p1.upper),]


  }

  for(i in 1:shells){
    left[[i]]<-peels[[i]][[1]]
    right[[i]]<-peels[[i]][[2]]

    df1 <- do.call("rbind", left)
    df2 <- do.call("rbind", right)
  }

  if(plot==TRUE){
    points(df1$x,df1$y,col="red", pch=16)
    points(df2$x,df2$y,col="blue", pch=16)}

  pointL<-df1
  pointR<-df2

  ## Clustering metrics-------------------------------------------------------------------

  #### 1. Calculating the euclidean distance of vertices to center------------------------

  ED1_sd<-sd(sqrt((mean(x)-df1[,1])^2+(mean(y)-df1[,2])^2))
  ED2_sd<-sd(sqrt((mean(x)-df2[,1])^2+(mean(y)-df2[,2])^2))

  #### 2. Calculating the perimeter and of the boundary points----------------------------

  perimL <- AP(concaveman(as.matrix(df1), concavity = alpha, length_threshold = 0))$Perimeter # left section
  perimR <-AP(concaveman(as.matrix(df2), concavity = alpha, length_threshold = 0))$Perimeter# right section

  areaL <- AP(concaveman(as.matrix(df1), concavity = alpha, length_threshold = 0))$Area # left section
  areaR <- AP(concaveman(as.matrix(df2), concavity = alpha, length_threshold = 0))$Area# right section


  ### Monte Carlo simulation for evidence testing  ---------------------------------------

  ED1_sim<-list()
  ED2_sim<-list()
  ED_all_sd_rise<-vector()
  ED_all_sd_fall<-vector()
  perimL_sim<-vector()
  perimR_sim<-vector()
  areaL_sim<-vector()
  areaR_sim<-vector()

  for(j in 1:simulations){

    peels<-list()
    left<-list()
    right<-list()

    ## simulation of data using summary statistics of the available data------------------

    set.seed(j)

    dat<-mvrnorm(n,mu=c(mean(x),mean(y)),Sigma)

    ## Determination of convex hull for the simulated data--------------------------------

    for(i in 1:shells){
      ch.index<-chull(dat) # find convex hull (index values for points)
      p1<-dat[ch.index,]   # extract peel 1
      p1_2<-p1[!duplicated(p1), ] # removes duplicate values i.e if two rows have same x and y values
      p1.upper<-upper.peel(p1_2)

      p1.upperx<-data.frame(x=p1.upper[,1],y=p1.upper[,2])

      ifelse( p1.upperx[1,2]!=max(p1.upperx$y),
              peels[[i]]<-split(p1.upperx,ifelse(p1.upperx$x<p1.upperx$x[which(p1.upperx$y==max(p1.upperx$y))],0,1)),
              peels[[i]]<-split(p1.upperx,ifelse(p1.upperx$x<= p1.upperx$x[which(p1.upperx$y==max(p1.upperx$y))],0,1)))


      index.func<-function(x,y){
        del<-list()
        for(i in 1:length(y[,1])){
          del[[i]]<-which(x[,1]==y[,1][i]& x[,2]==y[,2][i])
        }
        del

        d<-unlist(del, recursive = TRUE, use.names = TRUE) # convert all the indexed rows into one vector

        d
      }


      dat<-dat[-index.func(dat,p1.upper),]


    }

    for(i in 1:shells){
      left[[i]]<-peels[[i]][[1]]
      right[[i]]<-peels[[i]][[2]]

    }

    df1 <- do.call("rbind", left)
    df2 <- do.call("rbind", right)

    # 1. Enclidean distance measure-------------------------------------------------------

    ED1_2<-sqrt((mean(x)-df1[,1])^2+(mean(y)-df1[,2])^2)
    ED2_2<-sqrt((mean(x)-df2[,1])^2+(mean(y)-df2[,2])^2)

    ED1_sim[[j]]<-ED1_2
    ED2_sim[[j]]<-ED2_2

    # 2. Perimeter measure---------------------------------------------------------------

    perimL_2 <- AP(concaveman(as.matrix(df1), concavity = alpha, length_threshold = 0))$Perimeter # left section
    perimR_2 <- AP(concaveman(as.matrix(df2), concavity = alpha, length_threshold = 0))$Perimeter# right section
    areaL_2 <- AP(concaveman(as.matrix(df1), concavity = alpha, length_threshold = 0))$Area # left section
    areaR_2 <- AP(concaveman(as.matrix(df2), concavity = alpha, length_threshold = 0))$Area# right section

    perimL_sim[j]<-perimL_2
    perimR_sim[j]<-perimR_2
    areaL_sim[j]<-areaL_2
    areaR_sim[j]<-areaR_2
  }
  #---------------------------------------------------------------------------------------

  for(i in 1:simulations){
    ED_all_sd_rise[i]<-sd(ED1_sim[[i]])

  }



  for(i in 1:simulations){
    ED_all_sd_fall[i]<-sd(ED2_sim[[i]])
  }

  ## Calculating test indices-------------------------------------------------------------

  ### 1. sd test indices------------------------------------------------------------------

  p_sd_rise<-length(which(ED_all_sd_rise<=ED1_sd))/length(ED_all_sd_rise)
  p_sd_fall<-length(which(ED_all_sd_fall<=ED2_sd))/length(ED_all_sd_fall)

  MeanSDr<-mean(ED_all_sd_rise)
  MeanSDf<-mean(ED_all_sd_fall)

  ### 2. Perimeter area test indices-----------------------------------------------------------

  p_perim_rise<-pnorm(perimL, mean = mean(perimL_sim), sd = sd(perimL_sim))
  p_perim_fall<-pnorm(perimR, mean = mean(perimR_sim), sd = sd(perimR_sim))

  p_area_rise<-pnorm(areaL, mean = mean(areaL_sim), sd = sd(areaL_sim))
  p_area_fall<-pnorm(areaR, mean = mean(areaR_sim), sd = sd(areaR_sim))

  MeanperimL<-mean(perimL_sim)
  MeanperimR<-mean(perimR_sim)

  MeanareaL<-mean(areaL_sim)
  MeanareaR<-mean(areaR_sim)

  ## Output preparation-------------------------------------------------------------------

  if(method=="sd-enclidean"){

    Index<-c("sd","sd")
    Section<-c("Left","Right")
    value<-c(ED1_sd,ED2_sd)
    Mean<-c(MeanSDr,MeanSDf)
    p_value<-c(p_sd_rise, p_sd_fall)

    ## Plotting the data points for visualization-------------------------------------------

    if(plot==TRUE){
      hist(ED_all_sd_rise,freq = FALSE, xlab="sd",main = "Left")
      lines(density(ED_all_sd_rise), lwd = 1, col = "red")
      abline(v=ED1_sd, col="red",lty=2)

      hist(ED_all_sd_fall,freq = FALSE,xlab="sd",main = "Right")
      lines(density(ED_all_sd_fall), lwd = 1, col = "red")
      abline(v=ED2_sd, col="red",lty=2)
    }

  }

  if(method=="Perimeter"){

    Index<-c("Perimeter","Perimeter")
    Section<-c("Left","Right")
    value<-c(perimL,perimR)
    Mean<-c(MeanperimL,MeanperimR)
    p_value<-c(p_perim_rise, p_perim_fall)

    ## Plotting the data points for visualization-------------------------------------------

    if(plot==TRUE){

      hist(perimL_sim,freq = FALSE, xlab="Perimeter",main = "Left")
      lines(density(perimL_sim), lwd = 1, col = "red")
      abline(v=perimL, col="red",lty=2)

      hist(perimR_sim,freq = FALSE,xlab="Perimeter",main = "Right")
      lines(density(perimR_sim), lwd = 1, col = "red")
      abline(v=perimR, col="red",lty=2)
    }
  }

  if(method=="Area"){

    Index<-c("Area","Area")
    Section<-c("Left","Right")
    value<-c(areaL,areaR)
    Mean<-c(MeanareaL,MeanareaR)
    p_value<-c(p_area_rise, p_area_fall)

    ## Plotting the data points for visualization-----------------------------------------

    if(plot==TRUE){

      hist(areaL_sim,freq = FALSE, xlab="Area",main = "Left")
      lines(density(areaL_sim), lwd = 1, col = "red")
      abline(v=areaL, col="red",lty=2)

      hist(areaR_sim,freq = FALSE,xlab="Area",main = "Right")
      lines(density(areaR_sim), lwd = 1, col = "red")
      abline(v=areaR, col="red",lty=2)
    }
  }


  result<-data.frame(Index,Section,value, Mean,p_value)

  return(result)
}
