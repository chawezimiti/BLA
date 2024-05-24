#' Boundary line model determination using quantile regression
#'
#' This function fits a boundary model to the upper bounds of a scatter plot of
#' \code{x} and \code{y} by estimating the conditional quantile (0-1) of the
#' response variable, \code{y}, across values of the predictor variables, \code{x}.
#' This is achieved using optimization procedure and hence requires some starting
#' guess parameters of a proposed model.
#'
#' @param x A numeric vector of values for the independent variable.
#' @param y A numeric vector of values for the response variable.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for linear plateau model, \code{"mit"}
#'   for the Mitscherlich model, \code{"schmidt"} for the Schmidt model,
#'   \code{"logistic"} for logistic model, \code{"logisticND"} for logistic model
#'   proposed by Nelder (1961), \code{"inv-logistic"} for the inverse logistic
#'   model, \code{"double-logistic"} for the double logistic model, \code{"qd"}
#'   for quadratic model and the \code{"trapezium"} for the trapezium model.For custom
#'   models, set \code{model = "other"}.
#' @param equation A custom model function writen in the form of an R function. Applies
#'   only when argument \code{model="other"}, else it is \code{NULL}.
#' @param tau The quantile value (0- 1) that represents the boundary
#'   (\code{default is tau = 0.95}).
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
#' @param xmin Numeric value that describes the minimum \code{x} value to which the
#'   boundary line is to be fitted (default is \code{min(x)}).
#' @param xmax A numeric value that describes the maximum \code{x} value to which the
#'   boundary line is to be fitted (default is \code{max(x)}). \code{xmin} and
#'   \code{xmax} determine the subset of the data set used to fit boundary model.
#' @param line_smooth Parameter that describes the smoothness of the boundary line.
#'   (default is 1000). The higher the value, the smoother the line.
#' @param lwd Determines the thickness of the boundary line on the plot (default is 1).
#' @param line_col Selects the color of the boundary line.
#' @param optim.method Describes the method used to optimize the model as in the
#'   \code{optim()} function. The methods include \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} and \code{"Brent"}.
#' @param ... Additional graphical parameters.
#' @returns A list of length 5 consisting of the fitted model, equation form, parameters
#'   of the boundary line, the weighted residue sum square. Additionally, a graphical
#'   representation of the boundary line on the scatter plot is produced.
#'
#' @details
#' Some inbuilt models are available for the \code{blqr()} function. The suggest model
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
#'  \item Logistic model (\code{"logisticND"})  (Nelder (1961))
#'  \deqn{ y= \frac{\beta_0}{1+(\beta_1 \times e^{-\beta_2x})}}
#'   where \eqn{\beta_1} is a scaling parameter, \eqn{\beta_2} is a shape
#'   parameter and \eqn{\beta_0} is the maximum response.
#'
#'  \item Double logistic model (\code{"double-logistic"})
#'  \deqn{ y= \frac{\beta_{0,1}}{1+e^{\beta_2(\beta_1-x)}} -
#'  \frac{\beta_{0,2}}{1+e^{\beta_4(\beta_3-x)}}}
#'  where \eqn{\beta_1} is a scaling parameter one, \eqn{\beta_2} is a shape parameter
#'  one, \eqn{\beta_{0,1}} and \eqn{\beta_{0,2}} are the maximum response ,
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
#' settling for local optimum parameters, it is advised that the function is run
#' using several starting values and the results with the smallest error
#' (weighted residue sum square) can be taken as a representation of the global
#'  optimum.
#'
#'  The common errors encountered due to poor start values \enumerate{
#' \item function cannot be evaluated at initial parameters
#' \item initial value in 'vmmin' is not finite}
#'
#' @references
#'
#' Cade, B. S., & Noon, B. R. (2003). A gentle introduction to quantile regression
#' for ecologists. Frontiers in Ecology and the Environment, 1(8), 412-420.
#'
#' Nelder, J.A. 1961. The fitting of a generalization of the logistic curve.
#' Biometrics 17: 89–110.
#'
#' Phillips, B.F. & Campbell, N.A. 1968. A new method of fitting the von Bertelanffy
#' growth curve using data on the whelk. Dicathais, Growth 32: 317–329.
#'
#' Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line approach
#' to analyze N2O flux data from agricultural soils. Nutrient Cycling in Agroecosystems,
#' 57, 119-129.
#
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
#' @export
#'
#' @rdname blqr
#' @usage
#' blqr(x,y,model, equation=NULL,start,tau=0.95,optim.method="Nelder-Mead",
#'      xmin=min(bound$x),xmax=max(bound$x), plot=TRUE,line_col="red",lwd=1,
#'      line_smooth=1000,...)
#'
#' @examples
#'
#' x<-log(SoilP$P)
#' y<-SoilP$yield
#' start<-c(4,3,13.6)
#'
#' blqr(x,y, start=start,model = "lp", tau=0.99,
#'       xlab=expression("ET mm ha"^-1),
#'       ylab=expression("Wheat yield/ ton ha"^-1),
#'       pch=16, col="grey")
#'
blqr<-function(x,y,model, equation=NULL,start,tau=0.95,optim.method="Nelder-Mead",
               xmin=min(bound$x),xmax=max(bound$x),
               plot=TRUE,line_col="red",lwd=1,line_smooth=1000,...){

  #### Data preparation for quantile regression ------------------------------------------

  BLMod<-model
  if(plot==TRUE){plot(x,y,...)}

  bound <- na.omit(as.data.table(data.frame(x=x,y=y))) #removes NA's

  ## Setting data limits for boundary model fitting --------------------------------------

  L<-xmin
  U<-xmax

  if(L<min(bound$x)) stop("The set minimum limit is less than the mimum of bounding points")
  if(U>max(bound$x)) stop("The set maximum limit is greater than the maximum of bounding points")

  ifelse(L==min(bound$x), bound2<-bound, bound2<-bound[-which(bound$x<L),])
  ifelse(U==max(bound2$x), data1<-bound2, data1<-bound2[-which(bound2$x>U),])

  x<-data1$x
  y<-data1$y

  #### Fitting the two parameter Linear model --------------------------------------------

  if(model=="blm"){

    v<-length(start)
    if(v>2) stop("start has more than two values")
    if(v<2) stop("start has less than two values")

    trap<-function(x,ar,br){
      yr<-ar+br*x
      yout<-yr
      return(yout)
    }


    rss<-function(start,x,y){
      ar=start[1]
      br=start[2]

      yf<-unlist(lapply(x,FUN=trap,ar=ar,br=br))

      err<-(y-yf)
      errx<-sum(sum(err[which(err>0|err==0)])*tau +
                  sum(abs(err[which(err<0)]))*(1-tau))

      return(errx)
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

    ## Optimization using optim function--------------------------------------------------

    start=start[1:2]
    ooo<-optim(start,rss,x=x,y=y,hessian = T,method=optim.method)  #find LS estimate of start given data in x,yobs
    scale<-1/abs( parscale(ooo$par,x=x,y=y))
    oo<-optim(ooo$par,rss,x=x,y=y,hessian = T,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]

    if(plot==TRUE){
    xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth))
    yfit<-lapply(xfine,FUN=trap,ar=arf,br=brf)
    yfit<-unlist(yfit)

    lines(xfine,yfit,lwd=lwd,col=line_col)
    }

    hesmat<-oo$hessian
    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-oo$par
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082")

    RSS<-oo$value
    Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2082x")

    Parameters<-list(Model=BLMod,Equation=Equation,Parameters=estimates,RSS= RSS,Hessian= hesmat,
                       Start=start,optimMethod=optim.method,data=data1)

    class(Parameters) <- "cm" #necessary for only printing only part of the output

    return(Parameters)
  }

  #### Fitting the three parameter Linear model ------------------------------------------

  if(model=="lp"|model=="logistic"|model=="logisticND"|model=="inv-logistic"|model=="qd"|model=="mit"|model=="schmidt"){

    v<-length(start)
    if(v>3) stop("start has more than three values")
    if(v<3) stop("start has less than three values")

    ## set the function for each method---------------------------------------------------

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

    ## Loss function----------------------------------------------------------------------

    rss1<-function(start,x,y){
      ar=start[1]
      br=start[2]
      ym=start[3]

      yf<-unlist(lapply(x,FUN=trap1,ar=ar,br=br,ym=ym))

      err<-(y-yf)
      errx<-sum(sum(err[which(err>0|err==0)])*tau +
                  sum(abs(err[which(err<0)]))*(1-tau))

      return(errx)
    }

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

    ## Optimization using optim function--------------------------------------------------

    start=start[1:3]
    ooo<-optim(start,rss1,x=x,y=y, hessian = T,method=optim.method)
    scale<-1/abs( parscale1(ooo$par,x=x,y=y))
    oo<-optim(ooo$par,rss1,x=x,y=y, hessian = T,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo)

    arf=oo$par[1]
    brf=oo$par[2]
    ymf=oo$par[3]

    if(plot==TRUE){
      xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth))   #needs attention
      yfit<-lapply(xfine,FUN=trap1,ar=arf,br=brf,ym=ymf)
      yfit<-unlist(yfit)

      lines(xfine,yfit,lwd=lwd,col=line_col)}

    hesmat<-oo$hessian
    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-oo$par
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080")

    RSS<-oo$value

    Parameters<-structure(list(Model=BLMod,Equation=Equation,Parameters=estimates,RSS= RSS,Hessian= hesmat,
                     Start=start,optimMethod=optim.method,data=data1), class = "cm")

    class(Parameters) <- "cm" #necessary for only printing only part of the output
    return(Parameters)
  }

  #### Fitting the five parameter trapezium model ----------------------------------------

  if(model=="trapezium"){

    v<-length(start)
    if(v>5) stop("start has more than five values")
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

      err<-(y-yf)
      errx<-sum(sum(err[which(err>0|err==0)])*tau +
                  sum(abs(err[which(err<0)]))*(1-tau))

      return(errx)
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


    ## Optimization using the optim function----------------------------------------------

    start=start[1:5]
    ooo<-optim(start,rss2,x=x,y=y, hessian = T,method=optim.method)   #find LS estimate of start given data in x,yobs
    scale<-1/abs( parscale2(ooo$par,x=x,y=y))
    oo<-optim(ooo$par,rss2,x=x,y=y, hessian = T,method=optim.method,control = list(parscale = scale))

    ifelse(any(is.nan(oo$par))==T, oo<-ooo, oo<-oo) #rescalling sometimes produces NaN
    # and hence this make it to use the original values in ooo.

    arf=oo$par[1]
    brf=oo$par[2]
    ymf=oo$par[3]
    aff=oo$par[4]
    bff=oo$par[5]
    bp1<-(ymf-arf)/(brf) #estimates of the boundary break points
    bp2<-(ymf-aff)/(bff) #estimates of the boundary break points

    if(plot==TRUE){
      xfine=seq(min(x,na.rm = T),max(x,na.rm = T),(max(x,na.rm = T)-min(x,na.rm = T))/((max(x,na.rm = T)-min(x,na.rm = T))*line_smooth)) # this needs attention in all methods
      yfit<-lapply(xfine,FUN=trap2,ar=arf,br=brf,ym=ymf,af=aff,bf=bff)
      yfit<-unlist(yfit)

      lines(xfine,yfit,lwd=lwd,col=line_col)}

    hesmat<-oo$hessian
    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-oo$par
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080","\u03B2\u2083","\u03B2\u2084")


    RSS<-oo$value
    Equation<-noquote("y = min(\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080, \u03B2\u2083 + \u03B2\u2084x)")

    Parameters<-list(Model=BLMod,Equation=Equation,Parameters=estimates,RSS= RSS,Hessian= hesmat,
                     Start=start,optimMethod=optim.method,data=data1)

    class(Parameters) <- "cm" #necessary for only printing only part of the output
    return(Parameters)
  }

  #### Fitting the six parameter double-Logistic model -----------------------------------

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
      err<-(y-yf)
      errx<-sum(sum(err[which(err>0|err==0)])*tau +
                  sum(abs(err[which(err<0)]))*(1-tau))

      return(errx)
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

    ## Optimization using optim function--------------------------------------------------

    start=start[1:6]
    ooo<-optim(start,rss3,x=x,y=y, hessian = T,method=optim.method)
    scale<-1/abs( parscale3(ooo$par,x=x,y=y))
    oo<-optim(ooo$par,rss3,x=x,y=y, hessian = T,method=optim.method,control = list(parscale = scale))

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

      lines(xfine,yfit,lwd=lwd,col=line_col)
    }

    hesmat<-oo$hessian
    estimates<-matrix(NA,length(start),1,dimnames=list(c(),c("Estimate")))
    estimates[,1]<-c(arf,brf,ymf, ynf, aff, bff)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u20801","\u03B2\u20802","\u03B2\u2083","\u03B2\u2084")


    RSS<-oo$value
    Equation<-noquote("y = {\u03B2\u20801/1+[exp(\u03B2\u2082*(\u03B2\u2081-x))]} - {\u03B2\u20801/1+[exp(\u03B2\u2084*(\u03B2\u2083-x))]} ")

    Parameters<-list(Model=BLMod,Equation=Equation,Parameters=estimates,RSS= RSS,Hessian= hesmat,
                     Start=start,optimMethod=optim.method,data=data1)

    class(Parameters) <- "cm" #necessary for only printing only part of the output
    return(Parameters)
  }

  #### CUSTOM FUNCTIONS ------------------------------------------------------------------

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
      err<-(y-yf)
      errx<-sum(sum(err[which(err>0|err==0)])*tau +
                  sum(abs(err[which(err<0)]))*(1-tau))
      return(errx)
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

    ooo <- optim(start, rss4, x=x,y=y, method=optim.method, equation=equation)
    scale <- 1 / abs(parscale4(ooo$par, x=x,y=y, equation=equation))
    oo <- optim(ooo$par, rss4, x=x,y=y, method=optim.method, control=list(parscale=scale), equation=equation)

    if (any(is.nan(oo$par))) {
      oo <- ooo
    }

    param_values <- oo$par
    names(param_values) <- letters[1:length(param_values)]

    if (plot == TRUE) {
      xfine <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = line_smooth)
      yfit <- do.call(equation, c(list(x=xfine), as.list(param_values)))
      lines(xfine, yfit, lwd=lwd, col=line_col)
    }

    estimates <- matrix(param_values, nrow=length(param_values), ncol=1)
    rownames(estimates) <- names(param_values)
    colnames(estimates) <- "Estimate"

    RSS <- oo$value
    Equation<-equation # to print equation in output

    RSS<-oo$value

    Parameters<-structure(list(Model=BLMod,Equation=equation,Parameters=estimates,RSS= RSS), class = "cm")

    return(Parameters)

  }

}


