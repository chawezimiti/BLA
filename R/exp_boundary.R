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
#' @param plot If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
#'   is not part of output (default is \code{TRUE}).
#' @param ... Additional graphical parameters as with the \code{par()} function.
#'
#' @returns A dataframe with the p-values of obtaining the observed standard deviation
#'   of the euclidean distances of vertices in the upper peels to the center of the
#'   dataset for the left and right sections of the dataset.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
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
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' y<-evapotranspiration$`yield(t/ha)`
#' expl_boundary(x,y,10,1000)
#'
expl_boundary<-function(x,y,shells=10,simulations=1000,plot=TRUE,...){

  cat("Note: This function may take a few minutes to run for large datasets.\n\n")


  ## Selection of the x_min and x_max index values

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

  ## Removing NA'S from the data ######

  data<- data.frame(x=x,y=y)
  test<-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE)

  if(length(test)>0){
    data1<-data[-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE),]}else{
      data1<-data
    }

  x<-data1$x
  y<-data1$y
  dat<-cbind(x,y)

  ## setting the output area #######

  if(plot==TRUE){
    plot_layout<-rbind(c(1,1,2),c(1,1,3))
    layout(plot_layout)
    plot(dat,...)}

  ## Determination of the convex hull

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

  ## Calculating the euclidean distance of vertices to center

  ED1_sd<-sd(sqrt((mean(x)-df1[,1])^2+(mean(y)-df1[,2])^2))
  ED2_sd<-sd(sqrt((mean(x)-df2[,1])^2+(mean(y)-df2[,2])^2))


  ######### Monte Carlo simulation for evidence testing  #################################

  ED1_sim<-list()
  ED2_sim<-list()
  ED_all_sd_rise<-vector()
  ED_all_sd_fall<-vector()

  for(j in 1:simulations){

    peels<-list()
    left<-list()
    right<-list()

    # simulation of data using summary statistics of the available data

    dat<-mvrnorm(n,mu=c(mean(x),mean(y)),Sigma)
    x=dat[,1]
    y=dat[,2]

    ## Removal of outliers from the simulated data

    bag<-bagplot(x,y,na.rm = T,create.plot = FALSE)
    dat<-rbind(bag$pxy.bag,bag$pxy.outer)
    dat<-data.frame(x=dat[,1],y=dat[,2])

    ## Determination of convex hull for the simulated data

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


    ED1_2<-sqrt((mean(x)-df1[,1])^2+(mean(y)-df1[,2])^2)
    ED2_2<-sqrt((mean(x)-df2[,1])^2+(mean(y)-df2[,2])^2)

    ED1_sim[[j]]<-ED1_2
    ED2_sim[[j]]<-ED2_2

  }

  for(i in 1:simulations){
    ED_all_sd_rise[i]<-sd(ED1_sim[[i]])
  }



  for(i in 1:simulations){
    ED_all_sd_fall[i]<-sd(ED2_sim[[i]])
  }

  ## Calculating the sd test indices

  p_sd_rise<-length(which(ED_all_sd_rise<=ED1_sd))/length(ED_all_sd_rise)
  p_sd_fall<-length(which(ED_all_sd_fall<=ED2_sd))/length(ED_all_sd_fall)

  MeanSDr<-mean(ED_all_sd_rise)
  MeanSDf<-mean(ED_all_sd_fall)

  ## Output preparation

  Index<-c("sd","sd","Mean sd","Mean sd","p_value","p_value")
  value<-c(ED1_sd,ED2_sd,MeanSDr,MeanSDf,p_sd_rise, p_sd_fall)
  Section<-c("Left","Right","Left","Right","Left","Right")

  ## Plotting the data points for visualization

  if(plot==TRUE){
    hist(ED_all_sd_rise,freq = FALSE, xlab="sd",main = "Left")
    lines(density(ED_all_sd_rise), lwd = 1, col = "red")
    abline(v=ED1_sd, col="red",lty=2)

    hist(ED_all_sd_fall,freq = FALSE,xlab="sd",main = "Right")
    lines(density(ED_all_sd_fall), lwd = 1, col = "red")
    abline(v=ED2_sd, col="red",lty=2)

    par(mfrow=c(1,1))}

  data.frame(Index,Section,value)
}
