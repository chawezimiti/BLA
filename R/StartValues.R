#' Starting values for optimization functions
#'
#' This functions helps to determine initial values for a selected boundary line
#' model when using the functions \code{blbin()}, \code{blqr()}, \code{bolides()},
#' \code{cbvn()} and \code{ble_profile()} to determine model parameters.
#'
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for
#'   linear plateau model, \code{"mit"} for the Mitscherlich model, \code{"schmidt"}
#'   for the Schmidt model, \code{"logistic"} for logistic model, \code{"logisticND"}
#'   for logistic model proposed by Nelder (1961), \code{"inv-logistic"} for
#'   the inverse logistic model, \code{"double-logistic"} for the double logistic model,
#'   \code{"qd"} for quadratic model, \code{"trapezium"} for the trapezium model and
#'   \code{"explore"} for function use exploration. The default is \code{"explore"}.
#' @param p The number of selected points used to obtain start values for the logistic
#'   mitcherlich and schmidt models. It is \code{NULL} for other models.
#' @param digits Number of decimal points for logistic type models (default is 2).
#' @param ... Additional graphical parameters. Applies to the logistic, mitcherlich and
#'   schmidt models to control the text on the plot.
#' @details
#' This function uses the \code{locator()} function. Once the model is selected,
#' the points that make up the boundary points are selected using mouse click on
#' the plots.
#' @returns A list containing the parameters of the suggested model.
#'
#' @references
#'
#' Fekedulegn, D., Mac Siurtain, M.P., & Colbert, J.J. 1999. Parameter estimation of
#' nonlinear growth models in forestry. Silva Fennica 33(4): 327–336.
#'
#' Lark, R. M., & Milne, A. E. (2016). Boundary line analysis of the effect of water
#' filled pore space on nitrous oxide emission from cores of arable soil. European
#' Journal of Soil Science, 67 , 148-159.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#' @export
#' @examples
#' startValues(model="explore")
#'
#'
startValues<-function(model="explore",p=NULL,digits = 2,...){

  if (model != "blm" && model != "lp" && model != "qd" && model!="trapezium" &&
      model!="logistic" && model!="double-logistic" && model!="inv-logistic" &&
      model!="mit" && model!="schmidt" && model!="explore") {
    stop("model type not recorgnised ")
  }

  if(model=="blm"){
    n<-2
    cat("select 2 points on the plot that make up the linear model\n\n")
  }

  if(model=="lp"){
    n<-2
    cat("select 2 points on the plot that make up the linear-plateau model\n\n")
  }

  if(model=="qd"){
    n<-3
    cat("select 3 points on the plot that make up the quadratic model in ascending order of x\n\n")
  }

  if(model=="trapezium"){
    n<-4
    cat("select 4 points on the plot that make up the trapezium model in ascending order of x\n\n")
  }

  if(model=="logistic"){
    n<-p
    cat("select all points on the plot that make up the logistic model in ascending order of x\n\n")
  }

  if(model=="inv-logistic"){
    n<-p
    cat("select all points on the plot that make up the inv-logistic model in ascending order of x\n\n")
  }

  if(model=="double-logistic"){
    n<-p
    cat("select all points on the plot that make up the double-logistic model in ascending order of x\n\n")
  }

  if(model=="mit"){
    n<-p
    cat("select all points on the plot that make up the logistic model in ascending order of x\n\n")
  }

  if(model=="schmidt"){
    n<-p
    cat("select all points on the plot that make up the logistic model in ascending order of x\n\n")
  }


  if(model=="explore"){
    slope<-0
    intercept<-0
    cat("y = f(P\u2081,P\u2082|x)\n\n")

    names(intercept)<-c("P\u2081")
    names(slope)<-c("P\u2082")
    result<-list(Param1=intercept, Param2=slope)
    return(result)
  }

  if(model=="blm"){
    d<-locator(n)
    slope<-(d$y[2]-d$y[1])/(d$x[2]-d$x[1])
    intercept<- d$y[1]-slope*d$x[1]
    cat("y = \u03B2\u2081 + \u03B2\u2082*x\n\n")

    names(intercept)<-c("\u03B2\u2081")
    names(slope)<-c("\u03B2\u2082")
    result<-list(Intercept=intercept, slope=slope)
    return(result)
  }

  if(model=="lp"){
    d<-locator(n)
    slope<-(d$y[2]-d$y[1])/(d$x[2]-d$x[1])
    intercept<- d$y[1]-slope*d$x[1]
    maximum<-max(d$y[2],d$y[1])
    cat("y = min(\u03B2\u2081 + \u03B2\u2082*x, \u03B2\u2080)\n\n")

    names(intercept)<-c("\u03B2\u2081")
    names(slope)<-c("\u03B2\u2082")
    names(maximum)<-c("\u03B2\u2080")

    result<-list(max_response=maximum,Intercept=intercept, slope=slope)
    return(result)
  }

  if(model=="qd"){
    d<-locator(n)
    slope<-(d$y[2]-d$y[1])/(d$x[2]-d$x[1])
    intercept<- d$y[1]-slope*d$x[1]
    Co<-(d$y[2]-intercept-slope*d$x[2])/(d$x[2]*d$x[2])

    cat("\ny = \u03B2\u2081 + \u03B2\u2082x + \u03B2\u2083*x\u00B2\n\n")

    names(intercept)<-c("\u03B2\u2081")
    names(slope)<-c("\u03B2\u2082")
    names(Co)<-c("\u03B2\u2083")

    result<-list(Intercept=intercept, slope=slope, shape=Co)
    return(result)
  }

  if(model=="trapezium"){
    d<-locator(n)
    slope1<-(d$y[2]-d$y[1])/(d$x[2]-d$x[1])
    slope2<-(d$y[4]-d$y[3])/(d$x[4]-d$x[3])
    slopes=c(slope1,slope2)
    intercept1<- d$y[1]-slope1*d$x[1]
    intercept2<- d$y[3]-slope2*d$x[3]
    intercepts<-c(intercept1,intercept2)
    ymax<-mean(d$y[2],d$y[3])
    cat("y = min(\u03B2\u2081 + \u03B2\u2082*x, \u03B2\u2080, \u03B2\u2083 + \u03B2\u2084*x)\n\n")

    names(intercepts)<-c("\u03B2\u2081","\u03B2\u2083")
    names(slopes)<-c("\u03B2\u2082","\u03B2\u2084")
    names(ymax)<-c("\u03B2\u2080")

    result<-list(max_response=ymax,Intercepts=intercepts, slopes=slopes)
    return(result)
  }

  if(model=="logistic"){
    d<-locator(n)
    df<-data.frame(x=d$x,y=d$y)
    dx <- diff(df$x)
    dy <- diff(df$y)
    slopes <- dy / dx
    inflection<-median(df$x) # scaling value
    ymax<-max(df$y)

    x2<-vector()
    y2<-vector()
    for(i in 1:(length(df$x)-1)){
      x2[i]<-(df$x[i]+df$x[i+1])/2
      y2[i]<-(df$y[i]+df$y[i+1])/2
    }

    slopes2 <- round(slopes, digits = digits)
    text(x2,y2,slopes2, ...)

    cat("y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))\n\n")

    names(inflection)<-c("\u03B2\u2081")
    names(ymax)<-c("\u03B2\u2080")
    results<-list(max_response=ymax, scaling_parameter=inflection, shape_parameters =slopes2)
    return(results)

  }

  if(model=="inv-logistic"){
    d<-locator(n)
    df<-data.frame(x=d$x,y=d$y)
    dx <- diff(df$x)
    dy <- diff(df$y)
    slopes <- dy / dx
    inflection<-quantile(df$x, 0.75, names = FALSE) # scaling value
    ymax<-max(df$y)

    x2<-vector()
    y2<-vector()
    for(i in 1:(length(df$x)-1)){
      x2[i]<-(df$x[i]+df$x[i+1])/2
      y2[i]<-(df$y[i]+df$y[i+1])/2
    }

    slopes2 <- round(slopes, digits = digits)
    text(x2,y2,slopes2, ...)

    cat("y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))\n\n")

    names(inflection)<-c("\u03B2\u2081")
    names(ymax)<-c("\u03B2\u2080")
    results<-list(max_response=ymax, scaling_parameter=inflection, shape_parameters =slopes2)
    return(results)

  }

  if(model=="double-logistic"){
    d<-locator(n)
    df<-data.frame(x=d$x,y=d$y)
    dx <- diff(df$x)
    dy <- diff(df$y)
    slopes <- dy / dx
    inflection<-quantile(df$x, c(0.3,0.7), names = FALSE) # scaling value
    ymax<-max(df$y)
    max_response<-c(ymax,ymax)

    x2<-vector()
    y2<-vector()
    for(i in 1:(length(df$x)-1)){
      x2[i]<-(df$x[i]+df$x[i+1])/2
      y2[i]<-(df$y[i]+df$y[i+1])/2
    }

    slopes2 <- round(slopes, digits = digits)
    text(x2,y2,slopes2, ...)

    cat("y = \u03B2\u20801/(1+exp(\u03B2\u2082(\u03B2\u2081-x))) - \u03B2\u20802/(1+exp(\u03B2\u2084(\u03B2\u2083-x)))\n\n")

    names(inflection)<-c("\u03B2\u2081","\u03B2\u2083")
    names(max_response)<-c("\u03B2\u20801","\u03B2\u20802")

    results<-list(max_response=max_response , scaling_parameters=inflection, shape_parameters =slopes2)
    return(results)

  }

  if(model=="mit"){
    #Fekedulegn et al. (1999)

    d<-locator(n)
    df<-data.frame(x=d$x,y=d$y)
    ymax<-max(df$y) # Estimate of B_0

    dx <- diff(df$x)
    dy <- diff(df$y)
    slopes <- dy / dx

    shape_B2<- min(abs(slopes/ymax)) # Estimate of B_2 i.e rate of approach to B_0

    intercept<- df$y[2]-slopes[1]*df$x[2]

    B1<-(ymax-intercept)/min(abs(slopes)) # Estimate of magnitude of y at start of growth

    x2<-vector()
    y2<-vector()
    for(i in 1:(length(df$x)-1)){
      x2[i]<-(df$x[i]+df$x[i+1])/2
      y2[i]<-(df$y[i]+df$y[i+1])/2
    }

    slopes2 <- round(slopes/ymax, digits = digits)
    #text(x2,y2,slopes2, ...)

    cat("y = \u03B2\u2080 + \u03B2\u2081*\u03B2\u2082^x\n\n")

    names(B1)<-c("\u03B2\u2081")
    names(shape_B2)<-c("\u03B2\u2082")
    names(ymax)<-c("\u03B2\u2080")


    results<-list(max_response=ymax, scaling_parameter=B1, shape_parameter = shape_B2)
    return(results)

  }

  if(model=="schmidt"){
    #Lark, R. M., & Milne, A. E. (2016)
    d<-locator(n)
    df<-data.frame(x=d$x,y=d$y)
    ymax<-max(df$y) # Estimate of B_0

    dx <- diff(df$x)
    dy <- diff(df$y)
    slopes <- dy / dx

    shape_B2<- mean(df$x[which(df$y==max(df$y))]) # Estimate of B_2 i.e x at which y= maximum
                                                  # the mean is used in case there are two max values

    intercept<- df$y[2]-slopes[1]*df$x[2]

    B1<-(ymax-intercept)/min(abs(slopes))*-1 # Estimate of shape parameter. zero is a flat curve

    x2<-vector()
    y2<-vector()
    for(i in 1:(length(df$x)-1)){
      x2[i]<-(df$x[i]+df$x[i+1])/2
      y2[i]<-(df$y[i]+df$y[i+1])/2
    }

    slopes2 <- round(slopes/ymax, digits = digits)
    #text(x2,y2,slopes2, ...)

    cat("y = \u03B2\u2080 - \u03B2\u2081 (1-\u03B2\u2082)\u00B2)\n\n")

    names(B1)<-c("\u03B2\u2081")
    names(shape_B2)<-c("\u03B2\u2082")
    names(ymax)<-c("\u03B2\u2080")


    results<-list(max_response=ymax, scaling_parameter=B1,
                  shape_parameter = shape_B2)
    return(results)

  }


}
