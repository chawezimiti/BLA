#' Binning method for determining the boundary line model
#'
#' This function fits a boundary model to the upper bounds of a scatter plot of
#' \code{x} and \code{y} based on the binning method. The data are first divided
#' into equal sized sections in the x-axis and a boundary point in each section
#' is selected based on a set criteria (e.g. 0.90, 0.95 or 0.99 percentile of
#' \code{y} among other criteria). A model is then fitted to the resulting boundary
#' points by the least squares method. This is done using optimization procedure
#' and hence requires some starting guess parameters for the proposed model.
#'
#' @param x A numeric vector of values for the independent variable.
#' @param y A numeric vector of values for the response variable.
#' @param bins A numeric vector of length 3 or 4 that determines the size of sections.
#'   The first and second values give the range of the data to be binned while the
#'   third and fourth values give the width of the bins and the step size
#'   respectively. If only three values are provided, the step size is assumed to be
#'   equal to bin width.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"explore"} as default, \code{"blm"} for linear model, \code{"lp"} for
#'   linear plateau model, \code{"mit"} for the Mitscherlich model, \code{"schmidt"}
#'   for the Schmidt model, \code{"logistic"} for logistic model, \code{"logisticND"}
#'   for logistic model proposed by Nelder (1961), \code{"inv-logistic"} for
#'   the inverse logistic model, \code{"double-logistic"} for the double logistic model,
#'   \code{"qd"} for quadratic model and the \code{"trapezium"} for the trapezium model.
#'   The \code{"explore"} is used to check the position of boundary points in each bin
#'   so that the correct \code{model} can be applied. For custom models, set
#'   \code{model = "other"}.
#' @param equation A custom model function writen in the form of an R function. Applies
#'   only when argument \code{model="other"}, else it is \code{NULL}.
#'
#' @param start A numeric vector of initial starting values for optimization
#'   in fitting the boundary model. Its length and arrangement depend on the
#'   suggested model: \itemize{
#'   \item For the \code{"blm"} model, it is a vector of length 2 arranged as intercept
#'   and slope.
#'   \item For the \code{"lp"} model, it is a vector of length 3 arranged as intercept,
#'   slope and maximum response.
#'   \item For the \code{"logistic"} and \code{"inv-logistic"} models, it is a
#'   vector of length 3 arranged as the scaling parameter, shape parameter and maximum
#'   response.
#'   \item For the \code{"logisticND"} model proposed by Nelder (1961), it is a
#'   vector of length 3 arranged as the scaling parameter, shape parameter and maximum
#'   response.
#'   \item For the \code{"double-logistic"} model, it is a vector of length 6 arranged
#'   as the scaling parameter one, shape parameter one, maximum response, maximum
#'   response, scaling parameter two and shape parameter two.
#'   \item For the \code{"qd"} model, it is a vector of length 3 arranged as constant,
#'   linear coefficient and quadratic coefficient.
#'   \item For the \code{"trapezium"} model, it is a vector of length 3  arranged as
#'   intercept one, slope one, maximum response, intercept two and slope two.
#'   \item For the \code{"mit"} model, it is a vector of length 3 arranged as the
#'   intercept, shape parameter and the maximum response.
#'   \item For the \code{"schmidt"} model, it is a vector of length 3 arranged as scaling
#'   parameter, shape parameter (x-value at maximum response ) and maximum response.}
#' @param plot If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
#'   is not part of output (default is \code{TRUE}).
#' @param bp_pch Point character as \code{pch} of the \code{plot()} function. It controls
#'   the shape of the boundary points on plot (\code{bp_pch = 16} as default).
#' @param xmin Numeric value that describes the minimum \code{x} value to which the
#'   boundary line is to be fitted (default is \code{min(x)}).
#' @param xmax A numeric value that describes the maximum \code{x} value to which the
#'   boundary line is to be fitted (default is \code{max(x)}). \code{xmin} and
#'   \code{xmax} determine the subset of the data set used to fit boundary model.
#' @param line_smooth Parameter that describes the smoothness of the boundary line.
#'   (default is 1000). The higher the value, the smoother the line.
#' @param lwd Determines the thickness of the boundary line on the plot (default is 1).
#' @param bp_col Selects the color of the boundary points.
#' @param tau A percentile value (0-1) that represents the boundary point within each
#'   bin (default is \code{tau = 0.95}).
#' @param optim.method Describes the method used to optimize the model as in the
#'   \code{optim()} function. The methods include \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} and \code{"Brent"}.
#' @param bl_col Colour of the boundary line.
#' @param ... Additional graphical parameters as in the \code{par()} function.
#' @returns A list of length 5 consisting of the fitted model, equation form, parameters
#'   of the boundary line, the residue mean square and the boundary points. Additionally,
#'   a graphical representation of the boundary line on the scatter plot is produced.
#'
#' @details
#' Some inbuilt models are available for the \code{blbin()} function. The
#' \code{"explore"} option for the argument \code{model} generates a plot showing the
#' location of the boundary points selected by the binning procedure. This helps to
#' identify which model type is suitable to fit as a boundary line. The suggest model
#' forms are as follows: \enumerate{
#'  \item Linear model (\code{"blm"})
#'  \deqn{y=\beta_1 + \beta_2x}
#'  where \eqn{\beta_1} is the intercept and \eqn{\beta_2} is the slope.
#'
#'  \item Linear plateau model (\code{"lp"})
#'  \deqn{y= {\rm min}(\beta_1+\beta_2x, \beta_0)}
#'  where \eqn{\beta_1} is the intercept , \eqn{\beta_2} is the slope  and \eqn{\beta_0}
#'  is the maximum response.
#'
#'  \item The logistic (\code{"logistic"}) and inverse logistic (\code{"inv-logistic"})
#'  models
#'  \deqn{ y= \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}}
#'  \deqn{ y= \beta_0 - \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}}
#'  where \eqn{\beta_1} is a scaling parameter , \eqn{\beta_2} is a shape parameter
#'  and \eqn{\beta_0} is the maximum response.
#'
#'  \item Logistic model (\code{"logisticND"})  (Nelder (2009))
#'  \deqn{ y= \frac{\beta_0}{1+(\beta_1 \times e^{-\beta_2x})}}
#'   where \eqn{\beta_1} is a scaling parameter, \eqn{\beta_2} is a shape
#'   parameter and \eqn{\beta_0} is the maximum response.
#'
#'  \item Double logistic model (\code{"double-logistic"})
#'  \deqn{ y= \frac{\beta_{0,1}}{1+e^{\beta_2(\beta_1-x)}} -
#'  \frac{\beta_{0,2}}{1+e^{\beta_4(\beta_3-x)}}}
#'  where \eqn{\beta_1} is a scaling parameter one, \eqn{\beta_2} is shape parameter one,
#'  \eqn{\beta_{0,1}} and \eqn{\beta_{0,2}} are the maximum response ,
#'  \eqn{\beta_3} is a scaling parameter two and  \eqn{\beta_4} is a shape parameter two.
#'
#'  \item Quadratic model (\code{"qd"})
#'  \deqn{y=\beta_1 + \beta_2x + \beta_3x^2}
#'  where \eqn{\beta_1} is a constant, \eqn{\beta_2} is a linear coefficient
#'  and  \eqn{\beta_3} is the quadratic coefficient.
#'
#'  \item Trapezium model (\code{"trapezium"})
#'  \deqn{y={\rm min}(\beta_1+\beta_2x, \beta_0, \beta_3 + \beta_4x)}
#'  where  \eqn{\beta_1} is the intercept one, \eqn{\beta_2} is the slope one,
#'  \eqn{\beta_0} is the maximum response, \eqn{\beta_3} is the intercept two
#'  and \eqn{\beta_3} is the slope two.
#'
#'  \item Mitscherlich model (\code{"mit"})
#'  \deqn{y= \beta_0 - \beta_1*\beta_2^x}
#'  where \eqn{\beta_1} is the intercept, \eqn{\beta_2} is a shape parameter
#'  and \eqn{\beta_0} is the maximum response.
#'
#'  \item Schmidt model (\code{"schmidt"})
#'  \deqn{y= \beta_0 + \beta_1(x-\beta_2)^2}
#'  where \eqn{\beta_1} is ascaling parameter, \eqn{\beta_2} is a
#'  shape parameter (x-value at maximum response ) and \eqn{\beta_0} is the
#'  maximum response .
#'
#'  \item Custom model ("other")
#'  This option allows you to create your own model form using the function
#'  \code{function()}. The custom model should be assigned to the argument
#'  \code{equation}. Note that the parameters for the custom model should be
#'  \code{a} and \code{b} for a two parameter model; \code{a}, \code{b} and \code{c}
#'  for a three parameter model; \code{a}, \code{b}, \code{c} and \code{d} for a
#'  four parameter model and so on.
#'  }
#'
#' The function \code{blbin()} utilities the optimization procedure of the
#' \code{optim()} function to determine the model parameters. There is a tendency
#' for optimization algorithms to settle at a local optimum. To remove the risk of
#' settling for local optimum parameters, it is advised that the function is run using
#' several starting values and the results with the smallest error (residue mean square)
#' can be taken as a representation of the global optimum.
#'
#' The common errors encountered due to poor start values \enumerate{
#' \item function cannot be evaluated at initial parameters
#' \item initial value in 'vmmin' is not finite}
#'
#' @references
#' Casanova, D., Goudriaan, J., Bouma, J., & Epema, G. (1999). Yield gap analysis
#' in relation to soil properties in direct-seeded flooded rice.
#'
#' Nelder, J.A. 1961. The fitting of a generalization of the logistic curve.
#' Biometrics 17: 89–110.
#'
#' Phillips, B.F. & Campbell, N.A. 1968. A new method of fitting the von Bertelanffy
#' growth curve using data on the whelk. Dicathais, Growth 32: 317–329.
#'
#' Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line approach
#' to analyze N2O flux data from agricultural soils. Nutrient Cycling in Agro-ecosystems,
#' 57, 119-129.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
#' @export
#'
#' @rdname blbin
#' @usage
#' blbin(x,y,bins, model="explore", equation=NULL, start, tau=0.95,
#'       optim.method="Nelder-Mead", xmin=min(bound$x), xmax=max(bound$x),plot=TRUE,
#'       bp_col="red", bp_pch=16, bl_col="red", lwd=1,line_smooth=1000,...)
#'
#' @examples
#' x<-log(SoilP$P)
#' y<-SoilP$yield
#' start<-c(4,3,13.6, 35, -5)
#' bins<-c(1.6,4.74,0.314)
#'
#' blbin(x,y, bins=bins, start=start,model = "trapezium", tau=0.99,
#'        xlab=expression("Phosphorus/ln(mg L"^-1*")"),
#'        ylab=expression("Yield/ t ha"^-1), pch=16,
#'        col="grey", bp_col="grey")
#'
blbin<-function(x,y,bins,model="explore", equation=NULL,start, tau=0.95,
                optim.method="Nelder-Mead", xmin=min(bound$x),
                xmax=max(bound$x),plot=TRUE, bp_col="red", bp_pch=16, bl_col="red",
                lwd=1,line_smooth=1000, ...){

  BLMod<-model

########### DATA PREPARATION BINNING -----------------------------------------------------

##### Removing NA's ----------------------------------------------------------------------

  data<- data.frame(x=x,y=y)
  test<-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE)

  if(length(test)>0){
      df<-data[-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE),]}else{
      df<-data
    }

#### Determining the Bins ----------------------------------------------------------------

  if(length(bins) < 3) stop("The bins should be atleast length 3")
  if(length(bins) > 4) stop("The bins length should not be more than 4")


  if(length(bins)==3){     # For Fixed bins

    df<- df[order(df$x),]

    startpoint <- bins[1]  # start of window
    endpoint<-bins[2]      # end of window
    bin_size <- bins[3]    # window size
    step_size <- bins[3]   # Step size for sliding the window

    bound_y <- vector("numeric")  # Initialize vectors to store results
    avg_x <- vector("numeric")    # Initialize vectors to store results

    while (startpoint <= endpoint) {

      df1<-df[which(df$x > startpoint & df$x < startpoint+bin_size),]

      bound_y <- c(bound_y, as.numeric(quantile(df1$y, tau)))
      avg_x <- c(avg_x, mean(df1$x))
      startpoint <- startpoint + step_size
    }

    dataset<-data.frame(x=avg_x,y=bound_y)

    NANs<-c(which(is.nan( dataset$x)), which(is.nan( dataset$y))) # removing NaNs (non numbers)

    if(length(NANs)>0){
      dataset<- dataset[-NANs, ]
    }else{
      dataset <- dataset
    }
  }

  if(length(bins)==4){    # For sliding window

   df<- df[order(df$x),]

   startpoint <- bins[1]  # start point of window
   endpoint<-bins[2]      # start point of window
   bin_size <- bins[3]    # Window size
   step_size <- bins[4]   # Step size for sliding the window


   bound_y <- vector("numeric") # Initialize vectors to store results
   avg_x <- vector("numeric")   # Initialize vectors to store results

   while (startpoint <= endpoint) {

     df1<-df[which(df$x > startpoint & df$x < startpoint+bin_size),]

     bound_y <- c(bound_y, as.numeric(quantile(df1$y, tau)))
     avg_x <- c(avg_x, mean(df1$x))
     startpoint <- startpoint + step_size
   }

   dataset<-data.frame(x=avg_x,y=bound_y)

   NANs<-c(which(is.nan( dataset$x)), which(is.nan( dataset$y))) # removing NaNs (non numbers)

   if(length(NANs)>0){
     dataset<- dataset[-NANs, ]
   }else{
     dataset <- dataset
   }

 }

##### Plotting the boundary data for viewing ---------------------------------------------

  if(plot==TRUE){
    plot(x,y,...)
    points(dataset$x,dataset$y,col=bp_col, pch=bp_pch)
    }

#### Setting data limits for boundary fitting -------------------------------------------

  bound<-dataset

  L<-xmin
  U<-xmax

  if(L<min(bound$x)) stop("The set minimum limit is less than the mimum of bounding points")
  if(U>max(bound$x)) stop("The set maximum limit is greater than the maximum of bounding points")

  ifelse(L==min(bound$x), bound2<-bound, bound2<-bound[-which(bound$x<L),])
  ifelse(U==max(bound2$x), newdata5<-bound2, newdata5<-bound2[-which(bound2$x>U),])

#### Checking if boundary points have NA values ------------------------------------------

  test2<-which(is.na(newdata5$y)==TRUE)

  if(length(test2)>0) stop("Some bins do not contain response values. Make sure all bins contain data points")

#### Exploring data to choose model ------------------------------------------------------

  if(model=="explore"){
    if(plot==TRUE){
      plot(x,y,...)
      points(dataset$x,dataset$y, col=bp_col, pch=bp_pch)}
    return(summary(dataset))
  }

##### Fitting the two parameter linear model ---------------------------------------------

  if(model=="blm"){

    v<-length(start)
    if(v>2) warning("start has more than two values")
    if(v<2) stop("start has less than two values")

    trap<-function(x,ar,br){
      yr<-ar+br*x
      return(yr)
    }


    rss<-function(start,x,y){
      ar=start[1]
      br=start[2]

      yf<-unlist(lapply(x,FUN=trap,ar=ar,br=br))

      err<-sum((y-yf)^2)/length(x)
      return(err)
    }

    parscale<-function(a,x,y){

      eps=1e-4
      nr<-length(a)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(rss((a+del),x,y)-rss((a),x,y))/eps
      }

      return(part)
    }


    ## Optimization using optim function

    ooo<-optim(start,rss,x=newdata5$x,y=newdata5$y,method=optim.method)  #find LS estimate of start given data in x,yobs
    scale<-1/abs( parscale(ooo$par,x=newdata5$x,y=newdata5$y))
    oo<-optim(ooo$par,rss,x=newdata5$x,y=newdata5$y,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]

    if(plot==TRUE){ xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth))
    yfit<-lapply(xfine,FUN=trap,ar=arf,br=brf)
    yfit<-unlist(yfit)

    lines(xfine,yfit,lwd=lwd,col=bl_col)
    }

    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-c(arf,brf)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082")

    RMS<-oo$value
    Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2082x")

    Result<-list(Model=BLMod,Equation=Equation, Parameters=estimates, RMS=RMS, B.points=dataset)

    class(Result) <- "cm" #necessary for only printing only part of the output

    return(Result)
  }

#### Fitting the three parameter model ---------------------------------------------------

  if(model=="lp"|model=="logistic"|model=="logisticND"|model=="inv-logistic"|model=="qd"|model=="mit"|model=="schmidt"){

    v<-length(start)
    if(v>3) stop("start has more than three values")
    if(v<3) stop("start has less than three values")

    ## Set the function and loss function for each method

    if(model=="lp"){
      Equation<-noquote("y = min (\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080)")
      trap1<-function(x,ar,br,ym){
        yr<-ar+br*x
        yout<-min(c(yr,ym))
        return(yout)
      }}

    if(model=="logistic"){
      Equation<-noquote("y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))")
      trap1<-function(x,ar,br,ym){
        yr<-ym/(1+exp(br*(ar-x)))
        yout<-yr
        return(yout)
      }
    }

    if(model=="logisticND"){
      Equation<-noquote("y = \u03B2\u2080/1+[\u03B2\u2081exp(-\u03B2\u2082*x)]")
      trap1<-function(x,ar,br,ym){
        yr<-ym/(1+(ar*exp(-br*x)))
        yout<-yr
        return(yout)
      }

    }

    if(model=="inv-logistic"){
      Equation<-noquote("y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))")
      trap1<-function(x,ar,br,ym){
        yr<-ym-(ym/(1+exp(br*(ar-x))))
        yout<-yr
        return(yout)
      }
    }

    if(model=="qd"){
      Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2082x + \u03B2\u2083x\u00B2")
      trap1<-function(x,ar,br,ym){
        yr<- ar + br*x + ym*x*x
        yout<-yr

        return(yout)
      }
    }

    if(model=="mit"){
      Equation<-noquote("y = \u03B2\u2080 + \u03B2\u2081*\u03B2\u2082^x")
      trap1<-function(x,ar,br,ym){
        yr<-ym-ar*br^x
        yout<-yr
        return(yout)
      }
    }

    if(model=="schmidt"){
      Equation<-noquote("y = \u03B2\u2080 - \u03B2\u2081 (1-\u03B2\u2082)\u00B2)")
      trap1<-function(x,ar,br,ym){
        yr<-ym-ar*(x-br)*(x-br)
        yout<-yr
        return(yout)
      }
    }


    rss1<-function(start,x,y){
      ar=start[1]
      br=start[2]
      ym=start[3]

      yf<-unlist(lapply(x,FUN=trap1,ar=ar,br=br,ym=ym))

      err<-sum((y-yf)^2)/length(x)
      return(err)
    }

    ## scaling
    parscale1<-function(a,x,y){

      eps=1e-4
      nr<-length(a)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(rss1((a+del),x,y)-rss1((a),x,y))/eps
      }

      return(part)
    }


    ## Optimization using optim function

    ooo<-optim(start,rss1,x=newdata5$x,y=newdata5$y,method=optim.method)
    scale<-1/abs( parscale1(ooo$par,x=newdata5$x,y=newdata5$y))
    oo<-optim(ooo$par,rss1,x=newdata5$x,y=newdata5$y,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]
    ymf=oo$par[3]

    if(plot==TRUE){
      xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth))
      yfit<-lapply(xfine,FUN=trap1,ar=arf,br=brf,ym=ymf)
      yfit<-unlist(yfit)

      lines(xfine,yfit,lwd=lwd,col=bl_col)
    }


    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-c(arf,brf,ymf)
    if(model=="qd"){
      rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2083")}else{
        rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080")
      }

    RMS<-oo$value

    Parameters<-list(Model=BLMod,Equation=Equation,Parameters=estimates, RMS=RMS, Boundary_points= newdata5)
    class(Parameters) <- "cm" #necessary for only printing only part of the output

    return(Parameters)
  }

#### Fitting the five parameter Trapezium model ------------------------------------------

  if(model=="trapezium"){

    v<-length(start)
    if(v>5) warning("start has more than five values")
    if(v<5) stop("start has less than five values")


    trap2<-function(x,ar,br,ym,af,bf){
      yr<-ar+br*x
      yf<-af+bf*x

      yout<-min(c(yr,yf,ym))

      return(yout)
    }


    rss2<-function(start,x,y){
      ar=start[1]
      br=start[2]
      ym=start[3]
      af=start[4]
      bf=start[5]

      yf<-unlist(lapply(x,FUN=trap2,ar=ar,br=br,ym=ym,af=af,bf=bf))

      err<-sum((y-yf)^2)/length(x)
      return(err)
    }

    parscale2<-function(a,x,y){

      eps=1e-4
      nr<-length(a)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(rss2((a+del),x,y)-rss2((a),x,y))/eps
      }

      return(part)
    }


    ## Optimization

    ooo<-optim(start,rss2,x=newdata5$x,y=newdata5$y,method=optim.method)
    scale<-1/abs( parscale2(ooo$par,x=newdata5$x,y=newdata5$y))
    oo<-optim(ooo$par,rss2,x=newdata5$x,y=newdata5$y,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]
    ymf=oo$par[3]
    aff=oo$par[4]
    bff=oo$par[5]

    xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth)) # this needs attention in all methods
    yfit<-lapply(xfine,FUN=trap2,ar=arf,br=brf,ym=ymf,af=aff,bf=bff)
    yfit<-unlist(yfit)

    if(plot==TRUE){lines(xfine,yfit,lwd=lwd,col=bl_col)}

    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-c(arf,brf,ymf,aff,bff)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080","\u03B2\u2083","\u03B2\u2084")

    RMS<-oo$value
    Equation<-noquote("y = min(\u03B2\u2081+ \u03B2\u2082x, \u03B2\u2080, \u03B2\u2083 + \u03B2\u2084x)")

    Result<-list(Model=BLMod,Equation=Equation, Parameters=estimates, RMS=RMS, B.points=dataset)
    class(Result) <- "cm" #necessary for only printing only part of the output

    return(Result)
  }

#### Fitting the six parameter logistics model -------------------------------------------

  if(model=="double-logistic"){
    v<-length(start)

    if(v>6) warning("start has more than six values")
    if(v<6) stop("start has less than six values")


    trap3<-function(x,ar,br,ym,yn, af, bf){
      yr<-ym/(1 + exp((br*(ar-x)))) - yn/(1 + exp((bf*(af-x))))

      yout<-yr

      return(yout)
    }


    rss3<-function(start,x,y){
      ar=start[1]
      br=start[2]
      ym=start[3]
      yn=start[4]
      af=start[5]
      bf=start[6]

      yf<-unlist(lapply(x,FUN=trap3,ar=ar,br=br,ym=ym, yn=yn, af=af, bf=bf))

      err<-sum((y-yf)^2)/length(x)
      return(err)
    }

    parscale3<-function(a,x,y){

      eps=1e-4
      nr<-length(a)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(rss3((a+del),x,y)-rss3((a),x,y))/eps
      }

      return(part)
    }

    ## Optimization using optim function

    ooo<-optim(start,rss3,x=newdata5$x,y=newdata5$y,method=optim.method)
    scale<-1/abs( parscale3(ooo$par,x=newdata5$x,y=newdata5$y))
    oo<-optim(ooo$par,rss3,x=newdata5$x,y=newdata5$y,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]
    ymf=oo$par[3]
    ynf=oo$par[4]
    aff=oo$par[5]
    bff=oo$par[6]

    if(plot==TRUE){
      xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth))  #needs attention
      yfit<-lapply(xfine,FUN=trap3,ar=arf,br=brf,ym=ymf,yn=ynf,af=aff,bf=bff)
      yfit<-unlist(yfit)

      lines(xfine,yfit,lwd=lwd,col=bl_col)
    }

    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-c(arf,brf,ymf, ynf, aff, bff)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u20801","\u03B2\u20802","\u03B2\u2083","\u03B2\u2084")

    RMS<-oo$value
    Equation<-noquote("y = {\u03B2\u20801/1+[exp(\u03B2\u2082*(\u03B2\u2081-x))]} - {\u03B2\u20801/1+[exp(\u03B2\u2084*(\u03B2\u2083-x))]} ")

    Result<-list(Model=BLMod,Equation=Equation, Parameters=estimates, RMS=RMS, B.points=dataset)
    class(Result) <- "cm" #necessary for only printing only part of the output

    return(Result)
  }

#### USING CUSTOM FUNCTIONS --------------------------------------------------------------

  if(model=="other"){

    #### Names in start and rearranging them  --------------------------------------------

    are_entries_named <- function(vec) {
      # Check if names attribute is not NULL
      if (is.null(names(vec))) {
        return(FALSE)
      }

      # Check if all entries have non-NA and non-empty names
      has_valid_names <- all(!is.na(names(vec))) && all(names(vec) != "")
      return(has_valid_names)
    }

    if(are_entries_named(start)==TRUE){
      start<-start[order(names(start))]
    } else{
      start<-start
    }

    start<-unname(start) # removes names from start

    #### Dynamic parameter handling -------------------------------------------------------

    rss4 <- function(start, x, y, equation) {
      param_list <- as.list(start)
      names(param_list) <- letters[1:length(start)]
      yf <- do.call(equation, c(list(x=x), param_list))
      err <- sum((y - yf)^2) / length(x)
      return(err)
    }

    ## Scaling function for dynamic parameters --------------------------------------------

    parscale4 <- function(k, x, y, equation) {
      eps <- 1e-4
      nr <- length(k)
      part <- vector("numeric", nr)
      for (i in 1:nr) {
        del <- rep(0, nr)
        del[i] <- eps
        part[i] <- (rss4((k + del), x, y, equation) - rss4(k, x, y, equation)) / eps
      }
      return(part)
    }

    ## Optimization using optim function -------------------------------------------------

    ooo <- optim(start, rss4, x=newdata5$x, y=newdata5$y, method=optim.method, equation=equation)
    scale <- 1 / abs(parscale4(ooo$par, x=newdata5$x, y=newdata5$y, equation=equation))
    oo <- optim(ooo$par, rss4, x=newdata5$x, y=newdata5$y, method=optim.method, control=list(parscale=scale), equation=equation)

    if (any(is.nan(oo$par))) {
      oo <- ooo
    }

    param_values <- oo$par
    names(param_values) <- letters[1:length(param_values)]

    if (plot == TRUE) {
      xfine <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = line_smooth)
      yfit <- do.call(equation, c(list(x=xfine), as.list(param_values)))
      lines(xfine, yfit, lwd=lwd, col=bl_col)
    }

    estimates <- matrix(param_values, nrow=length(param_values), ncol=1)
    rownames(estimates) <- names(param_values)
    colnames(estimates) <- "Estimate"

    RMS <- oo$value
    Equation<-equation # to print equation in output

    Parameters <- list(Model=BLMod, Equation=equation, Parameters=estimates, RMS=RMS, Boundary_points=newdata5)
    class(Parameters) <- "cm" # necessary for only printing only part of the output

    return(Parameters)

  }

}
