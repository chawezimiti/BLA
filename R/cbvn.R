#' Fitting boundary line using censored bivariate normal model
#'
#' This function fits a response model to the upper limits  of a scatter plot of
#' of \code{x} and \code{y} to determine the most efficient response of \code{y}
#' as a function of \code{x} (given a measurement error of \code{y}) based on a
#' censored distribution (Milne et al., 2016). The location of censor in the data
#' cloud is determined based on the maximum likelihood approach. This is done using
#' optimization procedure and hence requires some starting guess parameters for the
#' proposed model. It then compares the results with an uncensored normal bivariate
#' distribution to access the appropriateness of the censored model.
#'
#' @param vals A dataframe with two numeric columns, independent (\code{x}) and
#'   dependent (\code{y}) variables respectively.
#' @param theta A numeric vector of initial starting values for optimization
#'   in fitting the boundary model. Its length and arrangement depend on the
#'   suggested model: \itemize{
#'   \item For the \code{"blm"} model, it is a vector of length 7 arranged as the
#'   intercept, the slope, mean of \code{x}, mean of \code{y}, standard deviation
#'   of \code{x}, standard deviation of \code{y} and the correlation of \code{x}
#'   and \code{y}.
#'   \item For the \code{"lp"} model, it is a vector of length 8 arranged as the
#'   intercept, the slope, the maximum or plateau response, mean of \code{x},
#'   mean of \code{y}, standard deviation of \code{x}, standard deviation of \code{y}
#'   and the correlation of \code{x} and \code{y}.
#'   \item For the \code{"mit"} model, it is a vector of length 8 arranged as the
#'   intercept, shape parameter, the maximum or plateau response, mean of \code{x},
#'   mean of \code{y}, standard deviation of \code{x}, standard deviation of \code{y}
#'   and the correlation of \code{x} and \code{y}.
#'   \item For the \code{"logistic"}, inv-logistic and "logisticfm" models,
#'   it is a vector of length 8 arranged as scaling parameter, shape parameter,
#'   the maximum or plateau value, mean of \code{x}, mean of \code{y},
#'   standard deviation of \code{x}, standard deviation of \code{y} and the
#'   correlation of \code{x} and \code{y}.
#'   \item For the \code{"double-logistic"} model,it is a vector of length 11
#'   arranged as scaling parameter, shape parameter, maximum response, maximum response,
#'   scaling parameter two, shape parameter two, mean of \code{x}, mean of \code{y},
#'   standard deviation of \code{x}, standard deviation of \code{y} and the correlation
#'   of \code{x} and \code{y}.
#'   \item For the \code{"trapezium"} model, it is a vector of length 10 arranged
#'   as intercept one, slope one, maximum response, intercept two, slope two,
#'   mean of \code{x}, mean of \code{y}, standard deviation of \code{x},
#'   standard deviation of \code{y} and the correlation of \code{x} and \code{y}.
#'   \item For the \code{"qd"} model, it is a vector of length 8 arranged as a
#'   constant, linear coefficient, quadratic coefficient,mean of \code{x},
#'   mean of \code{y}, standard deviation of \code{x}, standard deviation of \code{y}
#'   and the correlation of \code{x} and \code{y}.
#'   \item For the \code{"schmidt"} model, it is a vector of length 8 arranged the
#'   scaling parameter, shape parameter (x-value at maximum response ),
#'   maximum response, mean of \code{x}, mean of \code{y}, standard deviation of
#'   \code{x}, standard deviation of \code{y} and the correlation of \code{x}
#'   and \code{y}.
#' }
#' @param sigh Standard deviation of the measurement error.
#' @param UpLo Selects the type of boundary. \code{"U"} fits the upper boundary and
#'   "L" fits the lower boundary.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for linear plateau model, \code{"mit"}
#'   for the Mitscherlich model, \code{"schmidt"} for the Schmidt model, \code{"logistic"}
#'   for logistic model, \code{"logisticfm"} for logistic model proposed by
#'   Fermont et al (2009), \code{"inv-logistic"} for the inverse logistic model,
#'   \code{"double-logistic"} for the double logistic model, \code{"qd"} for
#'   quadratic model and the \code{"trapezium"} for the trapezium model.
#' @param equation A custom model function writen in the form of an R function. Applies
#'   only when argument \code{model="other"}, else it is \code{NULL}.
#' @param optim.method Describes the method used to optimize the model as in the
#'   \code{optim()} function. The methods include \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} and \code{"Brent"}.
#' @param plot If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
#'   is not part of output (default is \code{TRUE}).
#' @param Hessian If \code{True}, the hessian matrix is part of the output
#'   (default is \code{FALSE}`).
#' @param line_smooth Parameter that describes the smoothness of the boundary line.
#'   (default is 1000). The higher the value, the smoother the line.
#' @param lwd Determines the thickness of the boundary line on the plot (default is 1).
#' @param l_col Selects the color of the boundary line.
#' @param ... Additional graphical parameters as in the \code{par()} function.
#' @returns A list of length 5 consisting of the fitted model, equation form,
#'   parameters of the boundary line, AIC (for boundary line model and a null model)
#'   and a hessian matrix.  Additionally, a graphical representation of the boundary
#'   line on the scatter plot is produced.
#'
#' @details
#' Some inbuilt models are available for the \code{cbvn()} function. The suggest model
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
#'  \item Logistic model (\code{"logisticfm"})  (Fermont et al. 2009)
#'  \deqn{ y= \frac{\beta_0}{1+(\beta_1 \times e^{-\beta_2x})}}
#'   where \eqn{\beta_1} is a scaling parameter, \eqn{\beta_2} is a shape
#'   parameter and \eqn{\beta_0} is the maximum response.
#'
#'  \item Double logistic model (\code{"double-logistic"})
#'  \deqn{ y= \frac{\beta_{0,1}}{1+e^{\beta_2(\beta_1-x)}} -
#'  \frac{\beta_{0,2}}{1+e^{\beta_4(\beta_3-x)}}}
#'  where \eqn{\beta_1} is a scaling parameter one, \eqn{\beta_2} is a shape parameter one,
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
#'  \deqn{y= \beta_0 + \beta_1(1-e^{\frac{-x}{\beta_2}})}
#'  where \eqn{\beta_1} is the intercept, \eqn{\beta_2} is a shape parameter
#'  and \eqn{\beta_0} is the maximum response.
#'
#'  \item Schmidt model (\code{"schmidt"})
#'  \deqn{y= \beta_0 + \beta_1(x-\beta_2)^2}
#'  where \eqn{\beta_1} is ascaling parameter, \eqn{\beta_2} is a
#'  shape parameter (x-value at maximum response ) and \eqn{\beta_0} is the
#'  maximum response .
#'  }
#'
#' The function \code{cbvn()} utilities the optimization procedure of the
#' \code{optim()} function to determine the model parameters. There is a tendency
#' for optimization algorithms to settle at a local optimum. To remove the risk of
#' settling for local optimum parameters, it is advised that the function is run using
#' several starting values and the results with the smallest likelihood (or AIC)
#' can be taken as a representation of the global optimum.
#'
#' @references
#'
#' Dhanoa, M. S., Sanderson, R., Cardenas, L. M., Shepherd, A., Chadwick, D. R.,
#' Powell, C. D., ... & France, J. (2022). Overview and application of the
#' Mitscherlich equation and its extensions to estimate the soil nitrogen pool
#' fraction associated with crop yield and nitrous oxide emission. Advances in
#' Agronomy, 174, 269-295.
#'
#' Fermont, A. M., Van Asten, P. J., Tittonell, P., Van Wijk, M. T., &
#' Giller, K. E. (2009). Closing the cassava yield gap: an analysis from smallholder
#' farms in East Africa. Field Crops Research, 112 (1), 24–36.
#'
#' Lark, R. M., & Milne, A. E. (2016). Boundary line analysis of the effect of water
#' filled pore space on nitrous oxide emission from cores of arable soil. European
#' Journal of Soil Science, 67 , 148-159.
#'
#' Lark, R. M., Gillingham, V., Langton, D., & Marchant, B. P. (2020). Boundary line
#' models for soil nutrient concentrations and wheat yield in national-scale datasets.
#' European Journal of Soil Science, 71 , 334-351.
#'
#' Milne, A. E., Ferguson, R. B., & Lark, R. M. (2006). Estimating a boundary line
#' model for a biological response by maximum likelihood.Annals of Applied Biology,
#' 149, 223–234.
#'
#' Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line approach
#' to analyze N2O flux data from agricultural soils. Nutrient Cycling in Agroecosystems,
#' 57, 119-129.
#'
#' @author \enumerate{
#' \item Chawezi Miti \email{chawezi.miti@@nottingham.ac.uk}
#' \item Richard Murray Lark \email{murray.lark@@nottingham.ac.uk}
#' }
#' @import mvtnorm
#' @export
#'
#' @rdname cbvn
#' @usage
#' cbvn(vals,model="lp", equation=NULL, theta, sigh, UpLo="U", optim.method="BFGS",Hessian=FALSE,
#'      plot=TRUE, line_smooth=100, lwd=2, l_col="red",...)
#'
#' @examples
#' x<-evapotranspiration$`ET(mm)`
#' y<-evapotranspiration$`yield(t/ha)`
#' vals<-data.frame(x,y)
#' guess<-c(0.5,0.02,mean(x),mean(y),sd(x),sd(y),cor(x,y))
#'
#' cbvn(vals,theta = guess, sigh = 0.4,model= "blm")
#'
cbvn<-function(vals, model="lp", equation=NULL, theta, sigh, UpLo="U", optim.method="BFGS",
               Hessian=FALSE, plot=TRUE, line_smooth=100, lwd=2, l_col="red",...){

  cat("Note: This function may take a few minutes to run for large datasets.\n\n")

  data<-data.frame(x=vals[,1],y=vals[,2])

  ###Removing NA's ####################

  test<-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE)

  if(length(test)>0){
    vals<-data[-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE),]}else{
    vals<-data
    }
  ########################################

  sigh<-sigh #measurement error
  UpLo<-UpLo # set UpLo to "U" when fitting an upper boundary and "L" for a lower.
  #BLModxx<-model


  if(model=="lp"|model=="mit"|model=="logistic"|model=="inv-logistic"|model=="schmidt"|model=="qd"|model=="logisticfm"){

    v<-length(theta)
    if(v>8) stop("theta has more than eight values")
    if(v<8) stop("theta has less than eight values")

    #########################################################################
    if(model=="lp"){
      lp<-function(x,beta0,beta1,beta2){
        return(min(beta0,beta1+beta2*x))
      }

      BLMod<-lp
    }


    if(model=="mit"){
      mit<-function(x,beta0,beta1,beta2){
        return(beta1+beta0*(1-exp(-x/beta2)))
      }

      BLMod<-mit
    }

    if(model=="logistic"){
      logistic<-function(x,beta0,beta1,beta2){
        return(beta0/(1+exp(beta2*(beta1-x))))
      }

      BLMod<-logistic
    }


    if(model=="inv-logistic"){
      inv_logistic<-function(x,beta0,beta1,beta2){
        return(beta0-(beta0/(1+exp(beta2*(beta1-x)))))
      }

      BLMod<-inv_logistic
    }

    if(model=="logisticfm"){
      logisticfm<-function(x,beta0,beta1,beta2){
        return(beta0/(1+beta1*exp(-x*beta2)))
      }

      BLMod<-logisticfm
    }

    if(model=="schmidt"){
      schmidt<-function(x,beta0,beta1,beta2){
        return(beta0-beta2*(x-beta1)*(x-beta1))
      }

      BLMod<-schmidt
    }

    if(model=="qd"){
      qd<-function(x,beta0,beta1,beta2){
        return(beta1+beta2*x+beta0*x*x)
      }

      BLMod<-qd
    }


    drawBL<-function(x,beta0,beta1,beta2,BLMod){

      #BLGen<-match.fun(BLMod)
      y<-sapply(x,BLMod,beta0=beta0,beta1=beta1,beta2=beta2)
      return(y)
    }



    ###################### LIKELIHOOD FUNCTIONS#########################
    ################################################################
    nll_mef<-function(pars,uplo,BLMod){
      #########################################################################
      # Returns nll for 3-paramter BL model, parameters in pars,uplo specifies
      # upper or lower boundary.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #
      # BLMod is set to determine the parametric form of the BL model
      #
      # BL_bs:  broken stick.  beta0 is maximum, beta1 is intercept,
      #	  beta2 is slope.
      #
      # BL_mit: Mitscherlich.  beta0 is intercept, beta1 is shape parameter,
      #	  beta2 is max response - beta0
      #
      # Three parameters as currently set up.
      beta0<-pars[3]
      beta1<-pars[1]
      beta2<-pars[2]
      mux<-pars[4]
      muy<-pars[5]
      sdx<-pars[6]
      sdy<-pars[7]
      rcorr<-pars[8]

      if(uplo=="U"){
        nliks<-apply(vals,1,jdensup,BLMod=BLMod,
                     beta0=beta0,beta1=beta1,beta2=beta2,sigh=sigh,mux=mux,muy=muy,
                     sdx=sdx,sdy=sdy,rcorr=rcorr)
        nll<--sum(nliks)
      }else{
        print("Error, not set up for lower boundary")
        stop
      }
      return(nll)
    }

    ###########################################################################

    par_nll_mef<-function(x,UpLo,BLMod){# rough partial derivative at x of nll

      eps=1e-4

      nr<-length(x)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_mef((x+del),UpLo,BLMod)-nll_mef((x),UpLo,BLMod))/eps
      }

      return(part)
    }
    #########################################################################

    jdensup<-function(X,BLMod,beta0,beta1,beta2,sigh,mux,muy,sdx,sdy,rcorr){

      # joint density of observed values x and y given beta0,beta1,beta2 as bl
      # parameters (plateau, intercept and slope of bounded linear model).
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      #BLGen<-match.fun(BLMod)

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-BLMod(x,beta0,beta1,beta2)

      fy_x<-coffcturb(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }

    #########################################################################

    coffcturb<-function(x,mu,sig,a,c,sigh){
      #########################################################################
      #
      # a is right censor, c is left censor.  Set either to Inf/-Inf
      #
      # Notation as in Turban webpage, except k is substituted for c


      k<-((mu-c)/sig)
      d<-((mu-a)/sig)
      alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
      beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
      gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(d)-pnorm(k)))

      com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
      com2<-gamma*exp(com1)
      f<-com2*(pnorm((x-a-alpha)/beta)-pnorm((x-c-alpha)/beta))

      #rescale for censored
      f<-f*pnorm(c,mu,sig)

      #add contribution at x from mass at c

      f<-f+(dnorm((x-c),0,sigh)*(1-pnorm(c,mu,sig)))

      return(f)
    }
    #########################################################################
    ###################################################################


    #########################################################################
    nllmvn<-function(pars){
      #########################################################################
      # all data in ur


      mux<-pars[1]
      muy<-pars[2]
      sdx<-pars[3]
      sdy<-pars[4]
      rcorr<-pars[5]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy

      Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
      nllmvn<-0

      lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
      nllmvn<--sum(lliks)

      return(nllmvn)
    }

    ########################################################################


    #########################################################################
    nll_mef_maxyield<-function(pars){
      #########################################################################
      #
      # Returns nll for simple flat upper BL model, parameter in pars.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #

      ymax<-pars[1]
      mux<-pars[2]
      muy<-pars[3]
      sdx<-pars[4]
      sdy<-pars[5]
      rcorr<-pars[6]

      nliks<-apply(vals,1,jdens_maxyield,ymax=ymax,
                   sigh=sigh,mux=mux,muy=muy,
                   sdx=sdx,sdy=sdy,rcorr=rcorr)

      nll<--sum(nliks)

      return(nll)
    }


    #########################################################################
    jdens_maxyield<-function(X,ymax,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given a fixed maximum y
      # as the only model parameter.
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-ymax

      fy_x<-coffcturb(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }

    ####################END OF SUPPORT FUNCTIONS###################

    ############ OPTIMISING THE MODEL###############################

    mlest<-suppressWarnings(optim(theta,nll_mef,uplo=UpLo,BLMod=BLMod,
                                  method=optim.method,
                                  hessian="T"))

    scale<-suppressWarnings(1/abs(par_nll_mef(mlest$par,UpLo,BLMod)))

    mlest2<-suppressWarnings(optim(mlest$par,nll_mef,uplo=UpLo,BLMod=BLMod,
                                   method=optim.method,
                                   control = list(parscale = scale),
                                   hessian="T"))

    # there will be warning messages on running this as the optimizer will drift
    # into areas where density is very small.

    # compute the Akaike Information Criterion for the fitted model
    # This is 2* number of parameters + 2*(negative log likelihood)
    # number of parameters is 8

    AICbl<-(2*8)+(2*mlest2$value)

    # extract parameters of the boundary line

    beta0<-mlest2$par[3]
    beta1<-mlest2$par[1]
    beta2<-mlest2$par[2]

    ##Plotting the data

    if(plot==TRUE){ plot(vals,...)
      # set up values of the predictor variable to draw the boundary line
      #xdraw<-seq(1.5,4.5,0.01)
      x<-vals[,1]
      y<-vals[,2]
      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL((xdraw),beta0,beta1,beta2,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)
    }

    # Now compute the standard error of the boundary line model parameters then
    # write these out with the estimates

    hesmat<-mlest2$hessian

    estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
    estimates[,1]<-mlest2$par
    estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)

    if(model=="qd"){rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2083","mux","muy","sdx","sdy","rcorr")}else{

      rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080","mux","muy","sdx","sdy","rcorr")
    }

    #  Now fit the null model, an unbounded multivariate normal, and compute its
    # Akaike information criterion.

    theta2<-c(theta[c(4,5,6,7)],0)

    mvnmlest<-suppressWarnings(optim( theta2,nllmvn, #guess2 was c(1.6,9,0.5,1.7,0.0)
                                      method=optim.method))

    AImvn<-(2*5)+(2*mvnmlest$value) #


    #  Now fit a model with a constant upper boundary

    ymbest<-mlest2$par[3]+mlest2$par[2]
    ymax_guess<-c(ymbest,mlest2$par[c(4,5,6,7,8)])

    max_yield<-suppressWarnings(optim(c(ymax_guess),nll_mef_maxyield,
                                      method=optim.method))

    AImax_yield<-(2*6)+(2*max_yield$value)#

    #print(c(AImax_yield, AImvn,AICbl))

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")


    if(model=="lp"){ Equation<-noquote("y = min (\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080)")}
    if(model=="mit"){Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2080(1-exp(-x/\u03B2\u2081))")}
    if(model=="BL_logistic"){ Equation<-noquote("y = \u03B2\u2080/1+[\u03B2\u2081exp(-\u03B2\u2082x)]")}
    if(model=="BL_inv_logistic"){ Equation<-noquote("y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))")}
    if(model=="BL_logisticfm"){ Equation<-noquote("y = \u03B2\u2080/1+[\u03B2\u2081exp(-\u03B2\u2082*x)]")}
    if(model=="schmidt"){Equation<-noquote("y = \u03B2\u2080 - \u03B2\u2081 (1-\u03B2\u2082)\u00B2)")}
    if(model=="qd"){Equation<-noquote("y = \u03B2\u2081+\u03B2\u2082x+\u03B2\u2083x\u00B2")}

    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
    class(result)<-"cm"
    return(result)


  }

  ################### LINEAR MODEL ########################################################################################

  if(model=="blm"){

    v<-length(theta)
    if(v>7) stop("theta has more than seven values")
    if(v<7) stop("theta has less than seven values")

    ## The linear model defined ##

    blm<-function(x,beta0,beta1){
      return(beta0+beta1*x)
    }

    BLMod <- blm


    drawBL2<-function(x,beta0,beta1,BLMod){
      y<-sapply(x,BLMod,beta0=beta0,beta1=beta1)
      return(y)
    }


    #############Support functions for linear model##########################
    #########################################################################
    nll_mef2<-function(pars,uplo,BLMod){
      #########################################################################
      # Returns nll for 3-paramter BL model, parameters in pars,uplo specifies
      # upper or lower boundary.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #
      # BLMod is set to determine the parametric form of the BL model
      beta0<-pars[1]
      beta1<-pars[2]
      mux<-pars[3]
      muy<-pars[4]
      sdx<-pars[5]
      sdy<-pars[6]
      rcorr<-pars[7]

      if(uplo=="U"){
        nliks<-apply(vals,1,jdensup2,BLMod=BLMod,
                     beta0=beta0,beta1=beta1,sigh=sigh,mux=mux,muy=muy,
                     sdx=sdx,sdy=sdy,rcorr=rcorr)
        nll<--sum(nliks)
      }else{
        print("Error, not set up for lower boundary")
        stop
      }
      return(nll)
    }

    ###########################################################################
    par_nll_mef2<-function(x,UpLo,BLMod){# rough partial derivative at x of nll
      ###########################################################################

      eps=1e-4

      nr<-length(x)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_mef2((x+del),UpLo,BLMod)-nll_mef2((x),UpLo,BLMod))/eps
      }

      return(part)
    }
    #########################################################################

    #########################################################################
    jdensup2<-function(X,BLMod,beta0,beta1,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given beta0,beta1,beta2 as bl
      # parameters (plateau, intercept and slope of bounded linear model).
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      #BLGen2<-match.fun(BLMod)

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-BLMod(x,beta0,beta1)

      fy_x<-coffcturb(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }

    #########################################################################
    coffcturb<-function(x,mu,sig,a,c,sigh){
      #########################################################################
      #
      # a is right censor, c is left censor.  Set either to Inf/-Inf
      #
      # Notation as in Turban webpage, except k is substituted for c


      k<-((mu-c)/sig)
      d<-((mu-a)/sig)
      alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
      beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
      gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(d)-pnorm(k)))

      com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
      com2<-gamma*exp(com1)
      f<-com2*(pnorm((x-a-alpha)/beta)-pnorm((x-c-alpha)/beta))

      #rescale for censored
      f<-f*pnorm(c,mu,sig)

      #add contribution at x from mass at c

      f<-f+(dnorm((x-c),0,sigh)*(1-pnorm(c,mu,sig)))

      return(f)
    }
    #########################################################################
    ###################################################################


    #########################################################################
    nllmvn<-function(pars){
      #########################################################################
      # all data in ur


      mux<-pars[1]
      muy<-pars[2]
      sdx<-pars[3]
      sdy<-pars[4]
      rcorr<-pars[5]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy

      Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
      nllmvn<-0

      lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
      nllmvn<--sum(lliks)

      return(nllmvn)
    }

    ######################END OF SUPPORT FUNCTIONS##########################


    #### OPTIMISING THE LINEAR MODEL ########################

    mlest<-suppressWarnings(optim(theta,nll_mef2,uplo=UpLo,BLMod=BLMod,
                                  method=optim.method,
                                  hessian="T"))

    scale<-suppressWarnings(1/abs(par_nll_mef2(mlest$par,UpLo,BLMod)))

    mlest2<-suppressWarnings(optim(mlest$par,nll_mef2,uplo=UpLo,BLMod=BLMod,
                                   method=optim.method,
                                   control = list(parscale = scale),
                                   hessian="T"))

    # there will be warning messages on running this as the optimizer will drift
    # into areas where density is very small.

    # compute the Akaike Information Criterion for the fitted model
    # This is 2* number of parameters + 2*(negative log likelihood)
    # number of parameters is 8

    AICbl<-(2*8)+(2*mlest2$value)

    # extract parameters of the boundary line

    beta0<-mlest2$par[1]
    beta1<-mlest2$par[2]

    ##Plotting the data

    if(plot==TRUE){ plot(vals,...)

      x<-vals[,1]
      y<-vals[,2]

      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL2((xdraw),beta0,beta1,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)}

    hesmat<-mlest2$hessian
    #covpar<-solve(mlest2$hessian)


    estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
    estimates[,1]<-mlest2$par
    estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","mux","muy","sdx","sdy","rcorr")
    #print(estimates)

    #  Now fit the null model, an unbounded multivariate normal, and compute its
    # Akaike information criterion.

    theta2<-c(theta[c(3,4,5,6)],0)

    mvnmlest<- suppressWarnings(optim( theta2,nllmvn, #guess2 was c(1.6,9,0.5,1.7,0.0)
                                       method=optim.method))

    AImvn<-(2*5)+(2*mvnmlest$value) #


    #Output matrix
    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")

  if(model=="blm"){ Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2082x")}

  result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
  class(result)<-"cm"
  return(result)
  }


  ######################### TRAPEZIUM MODEL ###############################################################################

  if(model=="trapezium"){

    v<-length(theta)
    if(v>10) stop("theta has more than ten values")
    if(v<10) stop("theta has less than tenvalues")

    ### Define the trapezium function
    #########################################################################

    trapezium<-function(x,beta0,beta1,beta2,beta3,beta4){
      return(min(beta0,beta1+beta2*x,beta3+beta4*x))
    }

    BLMod <- trapezium

    drawBL3<-function(x,beta0,beta1,beta2,beta3,beta4,BLMod){
      y<-sapply(x,BLMod,beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
      return(y)
    }


    ######################LIKELIHOOD FUNCTIONS#########################
    ################################################################
    nll_mef3<-function(pars,uplo,BLMod){
      #########################################################################
      # Returns nll for 3-paramter BL model, parameters in pars,uplo specifies
      # upper or lower boundary.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #
      # BLMod is set to determine the parametric form of the BL model
      #
      # BL_bs:  broken stick.  beta0 is maximum, beta1 is intercept,
      #	  beta2 is slope.
      #
      # BL_mit: Mitscherlich.  beta0 is intercept, beta1 is shape parameter,
      #	  beta2 is max response - beta0
      #
      # Three parameters as currently set up.
      beta0<-pars[3]
      beta1<-pars[1]
      beta2<-pars[2]
      beta3<-pars[4]
      beta4<-pars[5]
      mux<-pars[6]
      muy<-pars[7]
      sdx<-pars[8]
      sdy<-pars[9]
      rcorr<-pars[10]

      if(uplo=="U"){
        nliks<-apply(vals,1,jdensup3,BLMod=BLMod,
                     beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,sigh=sigh,mux=mux,muy=muy,
                     sdx=sdx,sdy=sdy,rcorr=rcorr)
        nll<--sum(nliks)
      }else{
        print("Error, not set up for lower boundary")
        stop
      }
      return(nll)
    }

    ###########################################################################
    par_nll_mef3<-function(x,UpLo,BLMod){# rough partial derivative at x of nll
      ###########################################################################

      eps=1e-4

      nr<-length(x)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_mef3((x+del),UpLo,BLMod)-nll_mef3((x),UpLo,BLMod))/eps
      }

      return(part)
    }
    #########################################################################

    #########################################################################
    jdensup3<-function(X,BLMod,beta0,beta1,beta2,beta3,beta4,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given beta0,beta1,beta2 as bl
      # parameters (plateau, intercept and slope of bounded linear model).
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      #BLGen3<-match.fun(BLMod)

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-BLMod(x,beta0,beta1,beta2,beta3,beta4)

      fy_x<-coffcturb3(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }

    #########################################################################
    coffcturb3<-function(x,mu,sig,a,c,sigh){
      #########################################################################
      #
      # a is right censor, c is left censor.  Set either to Inf/-Inf
      #
      # Notation as in Turban webpage, except k is substituted for c


      k<-((mu-c)/sig)
      d<-((mu-a)/sig)
      alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
      beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
      gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(d)-pnorm(k)))

      com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
      com2<-gamma*exp(com1)
      f<-com2*(pnorm((x-a-alpha)/beta)-pnorm((x-c-alpha)/beta))

      #rescale for censored
      f<-f*pnorm(c,mu,sig)

      #add contribution at x from mass at c

      f<-f+(dnorm((x-c),0,sigh)*(1-pnorm(c,mu,sig)))

      return(f)
    }
    #########################################################################
    ###################################################################


    #########################################################################
    nllmvn3<-function(pars){
      #########################################################################
      # all data in ur


      mux<-pars[1]
      muy<-pars[2]
      sdx<-pars[3]
      sdy<-pars[4]
      rcorr<-pars[5]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy

      Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
      nllmvn<-0

      lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
      nllmvn<--sum(lliks)

      return(nllmvn)
    }

    ########################################################################


    #########################################################################
    nll_mef_maxyield3<-function(pars){
      #########################################################################
      #
      # Returns nll for simple flat upper BL model, parameter in pars.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #

      ymax<-pars[1]
      mux<-pars[2]
      muy<-pars[3]
      sdx<-pars[4]
      sdy<-pars[5]
      rcorr<-pars[6]

      nliks<-apply(vals,1,jdens_maxyield3,ymax=ymax,
                   sigh=sigh,mux=mux,muy=muy,
                   sdx=sdx,sdy=sdy,rcorr=rcorr)

      nll<--sum(nliks)

      return(nll)
    }


    #########################################################################
    jdens_maxyield3<-function(X,ymax,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given a fixed maximum y
      # as the only model parameter.
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-ymax

      fy_x<-coffcturb3(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }
    ############################

    ####################END OF SUPPORT FUNCTIONS###################
    ###############################################################

    ############OPTIMISING THE MODEL###############################

    mlest<-suppressWarnings(optim(theta,nll_mef3,uplo=UpLo,BLMod=BLMod,
                                  method=optim.method,
                                  hessian="T"))

    scale<-suppressWarnings(1/abs(par_nll_mef3(mlest$par,UpLo,BLMod)))

    mlest2<-suppressWarnings(optim(mlest$par,nll_mef3,uplo=UpLo,BLMod=BLMod,
                                   method=optim.method,
                                   control = list(parscale = scale),
                                   hessian="T"))

    # there will be warning messages on running this as the optimizer will drift
    # into areas where density is very small.

    # compute the Akaike Information Criterion for the fitted model
    # This is 2* number of parameters + 2*(negative log likelihood)
    # number of parameters is 8

    AICbl<-(2*8)+(2*mlest2$value)

    # extract parameters of the boundary line

    beta0<-mlest2$par[3]
    beta1<-mlest2$par[1]
    beta2<-mlest2$par[2]
    beta3<-mlest2$par[4]
    beta4<-mlest2$par[5]

    ##Plotting the data

    if(plot==TRUE){ plot(vals,...)

      x<-vals[,1]
      y<-vals[,2]
      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL3((xdraw),beta0,beta1,beta2,beta3,beta4,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)
    }

    # Now compute the standard error of the boundary line model parameters then
    # write these out with the estimates

    hesmat<-mlest2$hessian

    estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
    estimates[,1]<-mlest2$par

    estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080","\u03B2\u2083","\u03B2\u2084","mux","muy","sdx","sdy","rcorr")


    #  Now fit the null model, an unbounded multivariate normal, and compute its
    # Akaike information criterion.

    theta2<-c(theta[c(6,7,8,9)],0)

    mvnmlest<-suppressWarnings(optim( theta2,nllmvn3, #guess2 was c(1.6,9,0.5,1.7,0.0)
                                      method=optim.method))

    AImvn<-(2*5)+(2*mvnmlest$value) #


    #  Now fit a model with a constant upper boundary

    ymbest<-mlest2$par[3]+mlest2$par[2]
    ymax_guess<-c(ymbest,mlest2$par[c(6,7,8,9,10)])

    max_yield<-suppressWarnings(optim(c(ymax_guess),nll_mef_maxyield3,
                                      method=optim.method))

    AImax_yield<-(2*6)+(2*max_yield$value)#

    #print(c(AImax_yield, AImvn,AICbl))

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")


    if(model=="trapezium"){ Equation<-noquote("y = min (\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080, \u03B2\u2083 + \u03B2\u2084x)")}

    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
    class(result)<-"cm"
    return(result)


  }

  ##################  DOUBLE LOGISTIC MODEL ##############################################################################

  if(model=="double-logistic"){

    v<-length(theta)
    if(v>11) stop("theta has more than eleven values")
    if(v<11) stop("theta has less than eleven values")


    ###### Define the double logistic function ###########################

    double_logistic<-function(x,beta01,beta02,beta1,beta2,beta3,beta4){
      return((beta01/(1 + exp(beta2*(beta1-x)))) - (beta02/(1 + exp(beta4*(beta3-x)))))
    }

    BLMod <- double_logistic

    drawBL4<-function(x,beta01,beta02,beta1,beta2,beta3,beta4,BLMod){
      #BLGen4<-match.fun(BLMod)
      y<-sapply(x,BLMod,beta01=beta01,beta02=beta02,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
      return(y)
    }

    ######################LIKELIHOOD FUNCTIONS #########################

    nll_mef4<-function(pars,uplo,BLMod){

      beta1<-pars[1]
      beta2<-pars[2]
      beta01<-pars[3]
      beta02<-pars[4]
      beta3<-pars[5]
      beta4<-pars[6]
      mux<-pars[7]
      muy<-pars[8]
      sdx<-pars[9]
      sdy<-pars[10]
      rcorr<-pars[11]

      if(uplo=="U"){
        nliks<-apply(vals,1,jdensup4,BLMod=BLMod,
                     beta01=beta01,beta02=beta02,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,
                     sigh=sigh,mux=mux,muy=muy,
                     sdx=sdx,sdy=sdy,rcorr=rcorr)
        nll<--sum(nliks)
      }else{
        print("Error, not set up for lower boundary")
        stop
      }
      return(nll)
    }

    ###########################################################################

    par_nll_mef4<-function(x,UpLo,BLMod){# rough partial derivative at x of nll

      eps=1e-4

      nr<-length(x)
      part<-vector("numeric",nr)

      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_mef4((x+del),UpLo,BLMod)-nll_mef4((x),UpLo,BLMod))/eps
      }

      return(part)
    }
    #########################################################################

    #########################################################################
    jdensup4<-function(X,BLMod,beta01,beta02,beta1,beta2,beta3,beta4,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given beta0,beta1,beta2 as bl
      # parameters (plateau, intercept and slope of bounded linear model).
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      #BLGen4<-match.fun(BLMod)

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-BLMod(x,beta01,beta02,beta1,beta2,beta3,beta4)

      fy_x<-coffcturb4(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }
    #########################################################################

    coffcturb4<-function(x,mu,sig,a,c,sigh){

      # a is right censor, c is left censor.  Set either to Inf/-Inf
      #
      # Notation as in Turban webpage, except k is substituted for c


      k<-((mu-c)/sig)
      d<-((mu-a)/sig)
      alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
      beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
      gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(d)-pnorm(k)))

      com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
      com2<-gamma*exp(com1)
      f<-com2*(pnorm((x-a-alpha)/beta)-pnorm((x-c-alpha)/beta))

      #rescale for censored
      f<-f*pnorm(c,mu,sig)

      #add contribution at x from mass at c

      f<-f+(dnorm((x-c),0,sigh)*(1-pnorm(c,mu,sig)))

      return(f)
    }
    #########################################################################

    nllmvn4<-function(pars){
      #########################################################################
      # all data in ur


      mux<-pars[1]
      muy<-pars[2]
      sdx<-pars[3]
      sdy<-pars[4]
      rcorr<-pars[5]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy

      Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
      nllmvn<-0

      lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
      nllmvn<--sum(lliks)

      return(nllmvn)
    }

    ########################################################################

    nll_mef_maxyield4<-function(pars){
      #########################################################################
      #
      # Returns nll for simple flat upper BL model, parameter in pars.
      #
      # Measurement error fixed (sigh at top level)
      # Data in vals (n x 2 matrix) at top level.
      #

      ymax<-pars[1]
      mux<-pars[2]
      muy<-pars[3]
      sdx<-pars[4]
      sdy<-pars[5]
      rcorr<-pars[6]

      nliks<-apply(vals,1,jdens_maxyield4,ymax=ymax,
                   sigh=sigh,mux=mux,muy=muy,
                   sdx=sdx,sdy=sdy,rcorr=rcorr)

      nll<--sum(nliks)

      return(nll)
    }


    #########################################################################
    jdens_maxyield4<-function(X,ymax,sigh,mux,muy,sdx,sdy,rcorr){
      #########################################################################

      # joint density of observed values x and y given a fixed maximum y
      # as the only model parameter.
      # sigh ss measurement error
      # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
      # NB parameterization of rcorr to keep in [-1,1]

      x<-X[1]
      y<-X[2]

      rho<-tanh(rcorr)
      cov<-rho*sdx*sdy
      bet<-cov/(sdx*sdx)

      fx<-dnorm(x,mux,sdx)

      muyc<-muy+((x-mux)*bet)
      sdyc<-sdy*sqrt(1-(rho*rho))

      c<-ymax

      fy_x<-coffcturb4(y,muyc,sdyc,-Inf,c,sigh)

      fxy<-fy_x*fx
      return(log(fxy))
    }

    ####################END OF SUPPORT FUNCTIONS###################
    ###############################################################

    ############OPTIMISING THE MODEL###############################

    mlest<-suppressWarnings(optim(theta,nll_mef4,uplo=UpLo,BLMod=BLMod,
                                  method=optim.method,
                                  hessian="T"))

    scale<-suppressWarnings(1/abs(par_nll_mef4(mlest$par,UpLo,BLMod)))

    mlest2<-suppressWarnings(optim(mlest$par,nll_mef4,uplo=UpLo,BLMod=BLMod,
                                   method=optim.method,
                                   control = list(parscale = scale),
                                   hessian="T"))
    # there will be warning messages on running this as the optimizer will drift
    # into areas where density is very small.

    # compute the Akaike Information Criterion for the fitted model

    AICbl<-(2*8)+(2*mlest2$value)

    # Extract parameters of the boundary line
    beta1<-mlest2$par[1]
    beta2<-mlest2$par[2]
    beta01<-mlest2$par[3]
    beta02<-mlest2$par[4]
    beta3<-mlest2$par[5]
    beta4<-mlest2$par[6]

    # Plotting the data

    if(plot==TRUE){ plot(vals,...)
        x<-vals[,1]
        y<-vals[,2]
        xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
        ydraw<-drawBL4((xdraw),beta01,beta02,beta1,beta2,beta3,beta4,BLMod)
        lines(xdraw,ydraw,col=l_col,lwd=lwd)
    }

    # Now compute the standard error of the boundary line model parameters then
    # write these out with the estimates

    hesmat<-mlest2$hessian

    estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
    estimates[,1]<-mlest2$par

    estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u20801","\u03B2\u20802","\u03B2\u2083","\u03B2\u2084","mux","muy","sdx","sdy","rcorr")


    #  Now fit the null model, an unbounded multivariate normal, and compute its
    # Akaike information criterion.

    theta2<-c(theta[c(7,8,9,10)],0)

    mvnmlest<-suppressWarnings(optim( theta2,nllmvn4, #guess2 was c(1.6,9,0.5,1.7,0.0)
                                      method=optim.method))

    AImvn<-(2*5)+(2*mvnmlest$value) #


    #  Now fit a model with a constant upper boundary

    ymbest<-mlest2$par[3]+mlest2$par[1] # I changed here, it was mlest2$par[2]
    ymax_guess<-c(ymbest,mlest2$par[c(7,8,9,10,11)])

    max_yield<-suppressWarnings(optim(c(ymax_guess),nll_mef_maxyield4,
                                      method=optim.method))

    AImax_yield<-(2*6)+(2*max_yield$value)#

    #print(c(AImax_yield, AImvn,AICbl))

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")


    if(model=="BL_double_logistic"){
      Equation<-noquote("y = {\u03B2\u20801/1+[exp(\u03B2\u2082*(\u03B2\u2081-x))]} - {\u03B2\u20801/1+[exp(\u03B2\u2084*(\u03B2\u2083-x))]} ")}

    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
    class(result)<-"cm"
    return(result)


  }

  ######### CUSTOM MODELS ####################

  if(model=="other"){

      v<-length(theta)
      if(v>10) stop("Not set up for models with more than 5 parametrs. The argument theta should contain less than 10 values")
      Equation<-equation # to print equation in output
      theta<-unname(theta) # removes names from theta


      if(v==8){

        BLMod <- equation

        drawBL5 <- function(x,a,b,c,BLMod){
          y<-sapply(x,BLMod,a=a,b=b,c=c)
          return(y)
        }

        ###################### LIKELIHOOD FUNCTIONS#########################

        nll_mef5<-function(pars,uplo,BLMod){

          a<-pars[1]
          b<-pars[2]
          c<-pars[3]
          mux<-pars[4]
          muy<-pars[5]
          sdx<-pars[6]
          sdy<-pars[7]
          rcorr<-pars[8]

          if(uplo=="U"){
            nliks<-apply(vals,1,jdensup5,BLMod=BLMod,
                         a=a,b=b,c=c,sigh=sigh,mux=mux,muy=muy,
                         sdx=sdx,sdy=sdy,rcorr=rcorr)
            nll<--sum(nliks)
          }else{
            print("Error, not set up for lower boundary")
            stop
          }
          return(nll)
        }

        ###########################################################################

        par_nll_mef5<-function(x,UpLo,BLMod){# rough partial derivative at x of nll

          eps=1e-4

          nr<-length(x)
          part<-vector("numeric",nr)

          for (i in 1:nr){
            del<-rep(0,nr)
            del[i]<-eps
            part[i]<-(nll_mef5((x+del),UpLo,BLMod)-nll_mef5((x),UpLo,BLMod))/eps
          }

          return(part)
        }
        #########################################################################

        jdensup5<-function(X,BLMod,a,b,c,sigh,mux,muy,sdx,sdy,rcorr){

          # joint density of observed values x and y given beta0,beta1,beta2 as bl
          # parameters (plateau, intercept and slope of bounded linear model).
          # sigh ss measurement error
          # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
          # NB parameterization of rcorr to keep in [-1,1]

          x<-X[1]
          y<-X[2]

          rho<-tanh(rcorr)
          cov<-rho*sdx*sdy
          bet<-cov/(sdx*sdx)

          fx<-dnorm(x,mux,sdx)

          muyc<-muy+((x-mux)*bet)
          sdyc<-sdy*sqrt(1-(rho*rho))

          C<-BLMod(x,a,b,c)

          fy_x<-coffcturb5(y,muyc,sdyc,-Inf,C,sigh)

          fxy<-fy_x*fx
          return(log(fxy))
        }

        #########################################################################

        coffcturb5<-function(x,mu,sig,A,C,sigh){
          #########################################################################
          #
          # a is right censor, c is left censor.  Set either to Inf/-Inf
          #
          # Notation as in Turban webpage, except k is substituted for c


          k<-((mu-C)/sig)
          D<-((mu-A)/sig)
          alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
          beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
          gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(D)-pnorm(k)))

          com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
          com2<-gamma*exp(com1)
          f<-com2*(pnorm((x-A-alpha)/beta)-pnorm((x-C-alpha)/beta))

          #rescale for censored
          f<-f*pnorm(C,mu,sig)

          #add contribution at x from mass at c

          f<-f+(dnorm((x-C),0,sigh)*(1-pnorm(C,mu,sig)))

          return(f)
        }
        #########################################################################

        nllmvn5<-function(pars){

          mux<-pars[1]
          muy<-pars[2]
          sdx<-pars[3]
          sdy<-pars[4]
          rcorr<-pars[5]

          rho<-tanh(rcorr)
          cov<-rho*sdx*sdy

          Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
          nllmvn<-0

          lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
          nllmvn<--sum(lliks)

          return(nllmvn)
        }


        ####################END OF SUPPORT FUNCTIONS###################

        ############ OPTIMISING THE MODEL###############################

        mlest<-suppressWarnings(optim(theta,nll_mef5,uplo=UpLo,BLMod=BLMod,
                                      method=optim.method,
                                      hessian="T"))

        scale<-suppressWarnings(1/abs(par_nll_mef5(mlest$par,UpLo,BLMod)))

        mlest2<-suppressWarnings(optim(mlest$par,nll_mef5,uplo=UpLo,BLMod=BLMod,
                                       method=optim.method,
                                       control = list(parscale = scale),
                                       hessian="T"))

        # compute the Akaike Information Criterion for the fitted model


        AICbl<-(2*8)+(2*mlest2$value)

        # Extract parameters of the boundary line


        a<-mlest2$par[1]
        b<-mlest2$par[2]
        c<-mlest2$par[3]


        ## Plotting the data

        if(plot==TRUE){ plot(vals,...)
            x<-vals[,1]
            y<-vals[,2]
            xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
            ydraw<-drawBL5((xdraw),a,b,c,BLMod)
            lines(xdraw,ydraw,col=l_col,lwd=lwd)
        }

        # Now compute the standard error of the boundary line model parameters then
        # write these out with the estimates

        hesmat<-mlest2$hessian

        estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
        estimates[,1]<-mlest2$par
        estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)

        rownames(estimates)<-c("a","b","c","mux","muy","sdx","sdy","rcorr")

        #  Now fit the null model, an unbounded multivariate normal, and compute its
        # Akaike information criterion.

        theta2<-c(theta[c(4,5,6,7)],0)

        mvnmlest<-suppressWarnings(optim( theta2,nllmvn5, method=optim.method))

        AImvn<-(2*5)+(2*mvnmlest$value) #

        # Print the AIC values

        AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
        AikakeIC[,1]<-c(AImvn,AICbl)
        rownames(AikakeIC)<-c("mvn","BL")


        result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
        class(result)<-"cm"
        return(result)

      }

      if(v==9){

        BLMod <- equation

        drawBL6 <- function(x,a,b,c,d,BLMod){
          y<-sapply(x,BLMod,a=a,b=b,c=c,d=d)
          return(y)
        }

        ###################### LIKELIHOOD FUNCTIONS#########################

        nll_mef6<-function(pars,uplo,BLMod){

          a<-pars[1]
          b<-pars[2]
          c<-pars[3]
          d<-pars[4]
          mux<-pars[5]
          muy<-pars[6]
          sdx<-pars[7]
          sdy<-pars[8]
          rcorr<-pars[9]

          if(uplo=="U"){
            nliks<-apply(vals,1,jdensup6,BLMod=BLMod,
                         a=a,b=b,c=c,d=d,sigh=sigh,mux=mux,muy=muy,
                         sdx=sdx,sdy=sdy,rcorr=rcorr)
            nll<--sum(nliks)
          }else{
            print("Error, not set up for lower boundary")
            stop
          }
          return(nll)
        }

        ###########################################################################

        par_nll_mef6<-function(x,UpLo,BLMod){# rough partial derivative at x of nll

          eps=1e-4

          nr<-length(x)
          part<-vector("numeric",nr)

          for (i in 1:nr){
            del<-rep(0,nr)
            del[i]<-eps
            part[i]<-(nll_mef6((x+del),UpLo,BLMod)-nll_mef6((x),UpLo,BLMod))/eps
          }

          return(part)
        }
        #########################################################################

        jdensup6<-function(X,BLMod,a,b,c,d,sigh,mux,muy,sdx,sdy,rcorr){

          # joint density of observed values x and y given beta0,beta1,beta2 as bl
          # parameters (plateau, intercept and slope of bounded linear model).
          # sigh ss measurement error
          # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
          # NB parameterization of rcorr to keep in [-1,1]

          x<-X[1]
          y<-X[2]

          rho<-tanh(rcorr)
          cov<-rho*sdx*sdy
          bet<-cov/(sdx*sdx)

          fx<-dnorm(x,mux,sdx)

          muyc<-muy+((x-mux)*bet)
          sdyc<-sdy*sqrt(1-(rho*rho))

          C<-BLMod(x,a,b,c,d)

          fy_x<-coffcturb6(y,muyc,sdyc,-Inf,C,sigh)

          fxy<-fy_x*fx
          return(log(fxy))
        }

        #########################################################################

        coffcturb6<-function(x,mu,sig,A,C,sigh){
          #########################################################################
          #
          # a is right censor, c is left censor.  Set either to Inf/-Inf
          #
          # Notation as in Turban webpage, except k is substituted for c


          k<-((mu-C)/sig)
          D<-((mu-A)/sig)
          alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
          beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
          gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(D)-pnorm(k)))

          com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
          com2<-gamma*exp(com1)
          f<-com2*(pnorm((x-A-alpha)/beta)-pnorm((x-C-alpha)/beta))

          #rescale for censored
          f<-f*pnorm(C,mu,sig)

          #add contribution at x from mass at c

          f<-f+(dnorm((x-C),0,sigh)*(1-pnorm(C,mu,sig)))

          return(f)
        }
        #########################################################################

        nllmvn6<-function(pars){

          mux<-pars[1]
          muy<-pars[2]
          sdx<-pars[3]
          sdy<-pars[4]
          rcorr<-pars[5]

          rho<-tanh(rcorr)
          cov<-rho*sdx*sdy

          Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
          nllmvn<-0

          lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
          nllmvn<--sum(lliks)

          return(nllmvn)
        }


        ####################END OF SUPPORT FUNCTIONS###################

        ############ OPTIMISING THE MODEL###############################

        mlest<-suppressWarnings(optim(theta,nll_mef6,uplo=UpLo,BLMod=BLMod,
                                      method=optim.method,
                                      hessian="T"))

        scale<-suppressWarnings(1/abs(par_nll_mef6(mlest$par,UpLo,BLMod)))

        mlest2<-suppressWarnings(optim(mlest$par,nll_mef6,uplo=UpLo,BLMod=BLMod,
                                       method=optim.method,
                                       control = list(parscale = scale),
                                       hessian="T"))

        # compute the Akaike Information Criterion for the fitted model


        AICbl<-(2*8)+(2*mlest2$value)

        # Extract parameters of the boundary line


        a<-mlest2$par[1]
        b<-mlest2$par[2]
        c<-mlest2$par[3]
        d<-mlest2$par[4]


        ## Plotting the data

        if(plot==TRUE){ plot(vals,...)
          x<-vals[,1]
          y<-vals[,2]
          xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
          ydraw<-drawBL6((xdraw),a,b,c,d,BLMod)
          lines(xdraw,ydraw,col=l_col,lwd=lwd)
        }

        # Now compute the standard error of the boundary line model parameters then
        # write these out with the estimates

        hesmat<-mlest2$hessian

        estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
        estimates[,1]<-mlest2$par
        estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)

        rownames(estimates)<-c("a","b","c","d", "mux","muy","sdx","sdy","rcorr")

        #  Now fit the null model, an unbounded multivariate normal, and compute its
        # Akaike information criterion.

        theta2<-c(theta[c(4,5,6,7)],0)

        mvnmlest<-suppressWarnings(optim( theta2,nllmvn6, method=optim.method))

        AImvn<-(2*5)+(2*mvnmlest$value) #

        # Print the AIC values

        AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
        AikakeIC[,1]<-c(AImvn,AICbl)
        rownames(AikakeIC)<-c("mvn","BL")


        result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
        class(result)<-"cm"
        return(result)

      }

      if(v==10){

          BLMod <- equation
          drawBL7<-function(x,a,b,c,d,e,BLMod){
            y<-sapply(x,BLMod,a=a,b=b,c=c,d=d,e=e)
            return(y)
          }


          ######################LIKELIHOOD FUNCTIONS#########################
          ################################################################
          nll_mef7<-function(pars,uplo,BLMod){
            #########################################################################
            # Returns nll for 3-paramter BL model, parameters in pars,uplo specifies
            # upper or lower boundary.
            #
            # Measurement error fixed (sigh at top level)
            # Data in vals (n x 2 matrix) at top level.
            #
            # BLMod is set to determine the parametric form of the BL model
            #
            # BL_bs:  broken stick.  beta0 is maximum, beta1 is intercept,
            #	  beta2 is slope.
            #
            # BL_mit: Mitscherlich.  beta0 is intercept, beta1 is shape parameter,
            #	  beta2 is max response - beta0
            #
            # Three parameters as currently set up.

            a<-pars[1]
            b<-pars[2]
            c<-pars[3]
            d<-pars[4]
            e<-pars[5]
            mux<-pars[6]
            muy<-pars[7]
            sdx<-pars[8]
            sdy<-pars[9]
            rcorr<-pars[10]

            if(uplo=="U"){
              nliks<-apply(vals,1,jdensup7,BLMod=BLMod,
                           a=a,b=b,c=c,d=d,e=e,sigh=sigh,mux=mux,muy=muy,
                           sdx=sdx,sdy=sdy,rcorr=rcorr)
              nll<--sum(nliks)
            }else{
              print("Error, not set up for lower boundary")
              stop
            }
            return(nll)
          }

          ###########################################################################
          par_nll_mef7<-function(x,UpLo,BLMod){# rough partial derivative at x of nll
            ###########################################################################

            eps=1e-4

            nr<-length(x)
            part<-vector("numeric",nr)

            for (i in 1:nr){
              del<-rep(0,nr)
              del[i]<-eps
              part[i]<-(nll_mef7((x+del),UpLo,BLMod)-nll_mef7((x),UpLo,BLMod))/eps
            }

            return(part)
          }
          #########################################################################

          #########################################################################
          jdensup7<-function(X,BLMod,a,b,c,d,e,sigh,mux,muy,sdx,sdy,rcorr){
            #########################################################################

            # joint density of observed values x and y given beta0,beta1,beta2 as bl
            # parameters (plateau, intercept and slope of bounded linear model).
            # sigh ss measurement error
            # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
            # NB parameterization of rcorr to keep in [-1,1]

            x<-X[1]
            y<-X[2]

            #BLGen3<-match.fun(BLMod)

            rho<-tanh(rcorr)
            cov<-rho*sdx*sdy
            bet<-cov/(sdx*sdx)

            fx<-dnorm(x,mux,sdx)

            muyc<-muy+((x-mux)*bet)
            sdyc<-sdy*sqrt(1-(rho*rho))

            C<-BLMod(x,a,b,c,d,e)

            fy_x<-coffcturb7(y,muyc,sdyc,-Inf,C,sigh)

            fxy<-fy_x*fx
            return(log(fxy))
          }

          #########################################################################
          coffcturb7<-function(x,mu,sig,A,C,sigh){
            #########################################################################
            #
            # a is right censor, c is left censor.  Set either to Inf/-Inf
            #
            # Notation as in Turban webpage, except k is substituted for c


            k<-((mu-C)/sig)
            D<-((mu-A)/sig)
            alpha<-((sigh*sigh)*(x-mu))/((sigh*sigh)+(sig*sig))
            beta<-sqrt((sigh*sigh*sig*sig)/((sigh*sigh)+(sig*sig)))
            gamma<-(beta*sqrt(2*pi))/(2*pi*sig*sigh*(pnorm(D)-pnorm(k)))

            com1<--((x-mu)^2)/(2*((sigh*sigh)+(sig*sig)))
            com2<-gamma*exp(com1)
            f<-com2*(pnorm((x-A-alpha)/beta)-pnorm((x-C-alpha)/beta))

            #rescale for censored
            f<-f*pnorm(C,mu,sig)

            #add contribution at x from mass at c

            f<-f+(dnorm((x-C),0,sigh)*(1-pnorm(C,mu,sig)))

            return(f)
          }
          #########################################################################
          ###################################################################


          #########################################################################
          nllmvn7<-function(pars){
            #########################################################################
            # all data in ur


            mux<-pars[1]
            muy<-pars[2]
            sdx<-pars[3]
            sdy<-pars[4]
            rcorr<-pars[5]

            rho<-tanh(rcorr)
            cov<-rho*sdx*sdy

            Sigma<-matrix(c((sdx*sdx),cov,cov,(sdy*sdy)),2,2)
            nllmvn<-0

            lliks<-dmvnorm(vals, mean = c(mux,muy), sigma = Sigma, log = TRUE)
            nllmvn<--sum(lliks)

            return(nllmvn)
          }

          ########################################################################
#
#
#           #########################################################################
#           nll_mef_maxyield7<-function(pars){
#             #########################################################################
#             #
#             # Returns nll for simple flat upper BL model, parameter in pars.
#             #
#             # Measurement error fixed (sigh at top level)
#             # Data in vals (n x 2 matrix) at top level.
#             #
#
#             ymax<-pars[1]
#             mux<-pars[2]
#             muy<-pars[3]
#             sdx<-pars[4]
#             sdy<-pars[5]
#             rcorr<-pars[6]
#
#             nliks<-apply(vals,1,jdens_maxyield7,ymax=ymax,
#                          sigh=sigh,mux=mux,muy=muy,
#                          sdx=sdx,sdy=sdy,rcorr=rcorr)
#
#             nll<--sum(nliks)
#
#             return(nll)
#           }
#
#
#           #########################################################################
#           jdens_maxyield7<-function(X,ymax,sigh,mux,muy,sdx,sdy,rcorr){
#             #########################################################################
#
#             # joint density of observed values x and y given a fixed maximum y
#             # as the only model parameter.
#             # sigh ss measurement error
#             # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
#             # NB parameterization of rcorr to keep in [-1,1]
#
#             x<-X[1]
#             y<-X[2]
#
#             rho<-tanh(rcorr)
#             cov<-rho*sdx*sdy
#             bet<-cov/(sdx*sdx)
#
#             fx<-dnorm(x,mux,sdx)
#
#             muyc<-muy+((x-mux)*bet)
#             sdyc<-sdy*sqrt(1-(rho*rho))
#
#             c<-ymax
#
#             fy_x<-coffcturb7(y,muyc,sdyc,-Inf,c,sigh)
#
#             fxy<-fy_x*fx
#             return(log(fxy))
#           }
#           ############################

          ####################END OF SUPPORT FUNCTIONS###################
          ###############################################################

          ############OPTIMISING THE MODEL###############################

          mlest<-suppressWarnings(optim(theta,nll_mef7,uplo=UpLo,BLMod=BLMod,
                                        method=optim.method,
                                        hessian="T"))

          scale<-suppressWarnings(1/abs(par_nll_mef7(mlest$par,UpLo,BLMod)))

          mlest2<-suppressWarnings(optim(mlest$par,nll_mef7,uplo=UpLo,BLMod=BLMod,
                                         method=optim.method,
                                         control = list(parscale = scale),
                                         hessian="T"))

          # there will be warning messages on running this as the optimizer will drift
          # into areas where density is very small.

          # compute the Akaike Information Criterion for the fitted model
          # This is 2* number of parameters + 2*(negative log likelihood)
          # number of parameters is 8

          AICbl<-(2*8)+(2*mlest2$value)

          # extract parameters of the boundary line


          a<-mlest2$par[1]
          b<-mlest2$par[2]
          c<-mlest2$par[3]
          d<-mlest2$par[4]
          e<-mlest2$par[5]

          ##Plotting the data

          if(plot==TRUE){ plot(vals,...)

            x<-vals[,1]
            y<-vals[,2]
            xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
            ydraw<-drawBL7((xdraw),a,b,c,d,e,BLMod)
            lines(xdraw,ydraw,col=l_col,lwd=lwd)
          }

          # Now compute the standard error of the boundary line model parameters then
          # write these out with the estimates

          hesmat<-mlest2$hessian

          estimates<-matrix(NA,length(theta),2,dimnames=list(c(),c("Estimate","Standard error")))
          estimates[,1]<-mlest2$par

          estimates[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
          rownames(estimates)<-c("a","b","c","d","e","mux","muy","sdx","sdy","rcorr")


          #  Now fit the null model, an unbounded multivariate normal, and compute its
          # Akaike information criterion.

          theta2<-c(theta[c(6,7,8,9)],0)

          mvnmlest<-suppressWarnings(optim( theta2,nllmvn7, #guess2 was c(1.6,9,0.5,1.7,0.0)
                                            method=optim.method))

          AImvn<-(2*5)+(2*mvnmlest$value) #


          #  Now fit a model with a constant upper boundary
#
#           ymbest<-mlest2$par[3]+mlest2$par[2]
#           ymax_guess<-c(ymbest,mlest2$par[c(6,7,8,9,10)])
#
#           max_yield<-suppressWarnings(optim(c(ymax_guess),nll_mef_maxyield7,
#                                             method=optim.method))
#
#           AImax_yield<-(2*6)+(2*max_yield$value)#

          #print(c(AImax_yield, AImvn,AICbl))

          AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
          AikakeIC[,1]<-c(AImvn,AICbl)
          rownames(AikakeIC)<-c("mvn","BL")

          result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Hessian=hesmat)
          class(result)<-"cm"
          return(result)


      }

  }


}
