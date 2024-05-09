#' Likelihood profile for various measurement error values
#'
#' Estimates the standard deviation of measurement error (\code{sign}) of the response
#' variable, an input of the \code{cbvn()} function, when a measured value is not
#' available (Lark & Milne, 2016). \code{sigh} is fixed at each of a set of
#' values in turn, and remaining parameters are estimated conditional on \code{sigh}
#' by maximum likelihood. The maximized likelihoods for the sequence of values
#' constitutes a likelihood profile. The value of \code{sigh} where the profile
#' is maximized is selected.
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
#'   \item For the \code{"logistic"}, \code{"inv-logistic"} and \code{"logisticND"} models,
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
#' @param sigh A vector of the suggested standard deviations of the measurement error values.
#' @param UpLo Selects the type of boundary. \code{"U"} fits the upper boundary and
#'   "L" fits the lower boundary.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for linear plateau model, \code{"mit"}
#'   for the Mitscherlich model, \code{"schmidt"} for the Schmidt model, \code{"logistic"}
#'   for logistic model, \code{"logisticND"} for logistic model proposed by
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
#' @param ... Additional graphical parameters as in the \code{par()} function.
#' @returns  A list of length 2 containing the suggested standard deviations of
#'   measurement error values and the corresponding log-likelihood values.
#'   additionally, a likelihood profile plot (log-likelihood against the standard
#'   deviation of measurement error) is produced.
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
#'  \item Logistic model (\code{"logisticND"})  (Nelder (1961))
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
#' The function \code{ble_profile()} utilities the optimization procedure of the
#' \code{optim()} function to determine the model parameters. There is a tendency
#' for optimization algorithms to settle at a local optimum. To remove the risk of
#' settling for local optimum parameters, it is advised that the function is run
#' using several starting values and the results with the largest likelihood
#' can be taken as a representation of the global optimum.
#'
#' @references
#' Lark, R. M., & Milne, A. E. (2016). Boundary line analysis of the effect of water
#' filled pore space on nitrous oxide emission from cores of arable soil. European
#' Journal of Soil Science, 67 , 148-159.
#'
#' Nelder, J.A. 1961. The fitting of a generalization of the logistic curve.
#' Biometrics 17: 89–110.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
#' @export
#' @rdname ble_profile
#' @usage
#' ble_profile(vals, sigh, model="lp", equation=NULL,  theta, UpLo="U",
#'              optim.method="BFGS", plot=TRUE, ...)
#'
#' @examples
#'
#' x<-evapotranspiration$`ET(mm)`
#' y<-evapotranspiration$`yield(t/ha)`
#' vals<-data.frame(x,y)
#' theta<-c(0.5,0.02,289.6,2.39,83.8,1.05,0.295)
#' sigh<-c(0.08,0.1,0.11,0.13,0.15,0.2,0.3,0.4,0.5)
#' ble_profile(vals,sigh=sigh,theta=theta,model="blm")
#'
ble_profile<-function(vals, sigh, model="lp", equation=NULL, theta, UpLo="U", optim.method="BFGS", plot=TRUE, ...){

  cat("Note: This function may take a few minutes to run for large datasets.\n\n")

  data<-data.frame(x=vals[,1],y=vals[,2])
  ###Removing NA's
  test<-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE)

  if(length(test)>0){
    vals<-data[-which(is.na(data$x)==TRUE|is.na(data$y)==TRUE),]}else{
      vals<-data
    }
  #vals<-vals
  UpLo=UpLo
  BLMod<-model
  likelihood<-vector()


  if(model=="lp"|model=="mit"|model=="logistic"|model=="inv-logistic"|model=="logisticND"|model=="schmidt"|model=="qd"){

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

    if(model=="logisticND"){
      logisticND<-function(x,beta0,beta1,beta2){
        return(beta0/(1+beta1*exp(-x*beta2)))
      }

      BLMod<-logisticND
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



    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      ######################SUPPORT FUNCTIONS#########################
      ################################################################
      nll_mef<-function(pars,uplo,BLMod){

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
                       beta0=beta0,beta1=beta1,beta2=beta2,sigh=sigh[i],mux=mux,muy=muy,
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
        ###########################################################################

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

      #########################################################################
      jdensup<-function(X,BLMod,beta0,beta1,beta2,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
        #########################################################################

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

      coffcturb<-function(x,mu,sig,a,c,sigh=sigh[i]){

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


      ####################END OF SUPPORT FUNCTIONS###################
      ###############################################################

      ############OPTIMISING THE MODEL###############################
      #  Estimates are found by minimizing the negative log likelihood.  The first estimate
      #  starting from guess appears as mlest$par.  scale is recomputed from this, and
      #  a second run is started from this first solution

      mlest<-suppressWarnings(optim(theta,nll_mef,uplo=UpLo,BLMod=BLMod,
                                    method=optim.method,
                                    hessian="T"))

      scale<-suppressWarnings(1/abs(par_nll_mef(mlest$par,UpLo,BLMod)))

      mlest2<-suppressWarnings(optim(mlest$par,nll_mef,uplo=UpLo,BLMod=BLMod,
                                     method=optim.method,
                                     control = list(parscale = scale),
                                     hessian="T"))

      likelihood[i]<-mlest2$value*-1

    }


  }

  if(model=="blm"){

    v<-length(theta)
    if(v>7) stop("theta has more than seven values")
    if(v<7) stop("theta has less than seven values")


    blm<-function(x,beta0,beta1){
      return(beta0+beta1*x)
    }

    BLMod<-blm

    for(i in 1:length(sigh)){

      #sigh<-sigh[i]
      #likelihood<-vector()
      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      #############SUPPORT FUNCTIONS FOR LINEAR MODEL##########################
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
                       beta0=beta0,beta1=beta1,sigh=sigh[i],mux=mux,muy=muy,
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
      jdensup2<-function(X,BLMod,beta0,beta1,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
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
      coffcturb<-function(x,mu,sig,a,c,sigh=sigh[i]){
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


      ######################END OF SUPPORT FUNCTIONS##########################
      ########################################################################

      ####OPTIMISING THE LINEAR MODEL########################
      #######################################################
      #  Estimates are found by minimizing the negative log likelihood.  The first estimate
      #  starting from guess appears as mlest$par.  scale is recomputed from this, and
      #  a second run is started from this first solution

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


      likelihood[i]<-mlest2$value*-1

    }


  }

  if(model=="trapezium"){

    v<-length(theta)
    if(v>10) stop("theta has more than ten values")
    if(v<10) stop("theta has less than ten values")


    trapezium<-function(x,beta0,beta1,beta2,beta3,beta4){
      return(min(beta0,beta1+beta2*x,beta3+beta4*x))
    }

    BLMod <- trapezium

    for(i in 1:length(sigh)){

      #sigh<-sigh[i]
      #likelihood<-vector()
      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals
    ######################SUPPORT FUNCTIONS#########################
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
                     beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,sigh=sigh[i],mux=mux,muy=muy,
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
    jdensup3<-function(X,BLMod,beta0,beta1,beta2,beta3,beta4,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
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
    coffcturb3<-function(x,mu,sig,a,c,sigh=sigh[i]){
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
    ###############################################################


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

    likelihood[i]<-mlest2$value*-1
    }

  }

  if(model=="double-logistic"){

    v<-length(theta)
    if(v>11) stop("theta has more than eleven values")
    if(v<11) stop("theta has less than eleven values")


    double_logistic<-function(x,beta01,beta02,beta1,beta2,beta3,beta4){
      return((beta01/(1 + exp(beta2*(beta1-x)))) - (beta02/(1 + exp(beta4*(beta3-x)))))
    }

    BLMod<-double_logistic

    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      ######################SUPPORT FUNCTIONS#########################
      ################################################################
      nll_mef4<-function(pars,uplo,BLMod){
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
                       sigh=sigh[i],mux=mux,muy=muy,
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
        ###########################################################################

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
      jdensup4<-function(X,BLMod,beta01,beta02,beta1,beta2,beta3,beta4,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
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

      coffcturb4<-function(x,mu,sig,a,c,sigh=sigh[i]){
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

      likelihood[i]<-mlest2$value*-1
    }

  }

  if(model=="other"){

    v<-length(theta)
    if(v>10) stop("Not set up for models with more than 5 parametrs. The argument theta should contain less than 10 values")
    Equation<-equation # to print equation in output
    theta<-unname(theta) # removes names from theta

    if(v==8){

      BLMod <- equation

      for(i in 1:length(sigh)){

        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ######################SUPPORT FUNCTIONS#########################
        ################################################################
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
                         a=a,b=b,c=c,sigh=sigh[i],mux=mux,muy=muy,
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
          ###########################################################################

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

        #########################################################################
        jdensup5<-function(X,BLMod,a,b,c,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
          #########################################################################

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

          C<-BLMod(x,a,b,c)

          fy_x<-coffcturb5(y,muyc,sdyc,-Inf,C,sigh)

          fxy<-fy_x*fx
          return(log(fxy))
        }
        #########################################################################

        coffcturb5<-function(x,mu,sig,A,C,sigh=sigh[i]){

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


        ####################END OF SUPPORT FUNCTIONS###################
        ###############################################################

        ############OPTIMISING THE MODEL###############################
        #  Estimates are found by minimizing the negative log likelihood.  The first estimate
        #  starting from guess appears as mlest$par.  scale is recomputed from this, and
        #  a second run is started from this first solution

        mlest<-suppressWarnings(optim(theta,nll_mef5,uplo=UpLo,BLMod=BLMod,
                                      method=optim.method,
                                      hessian="T"))

        scale<-suppressWarnings(1/abs(par_nll_mef5(mlest$par,UpLo,BLMod)))

        mlest2<-suppressWarnings(optim(mlest$par,nll_mef5,uplo=UpLo,BLMod=BLMod,
                                       method=optim.method,
                                       control = list(parscale = scale),
                                       hessian="T"))

        likelihood[i]<-mlest2$value*-1

      }
    }


    if(v==9){

      BLMod <- equation

      for(i in 1:length(sigh)){

        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ######################SUPPORT FUNCTIONS#########################

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
                         a=a,b=b,c=c,d=d,sigh=sigh[i],mux=mux,muy=muy,
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
          ###########################################################################

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

        #########################################################################
        jdensup6<-function(X,BLMod,a,b,c,d,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
          #########################################################################

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

          C<-BLMod(x,a,b,c,d)

          fy_x<-coffcturb6(y,muyc,sdyc,-Inf,C,sigh)

          fxy<-fy_x*fx
          return(log(fxy))
        }
        #########################################################################

        coffcturb6<-function(x,mu,sig,A,C,sigh=sigh[i]){

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


        ####################END OF SUPPORT FUNCTIONS###################
        ###############################################################

        ############OPTIMISING THE MODEL###############################
        #  Estimates are found by minimizing the negative log likelihood.  The first estimate
        #  starting from guess appears as mlest$par.  scale is recomputed from this, and
        #  a second run is started from this first solution

        mlest<-suppressWarnings(optim(theta,nll_mef6,uplo=UpLo,BLMod=BLMod,
                                      method=optim.method,
                                      hessian="T"))

        scale<-suppressWarnings(1/abs(par_nll_mef6(mlest$par,UpLo,BLMod)))

        mlest2<-suppressWarnings(optim(mlest$par,nll_mef6,uplo=UpLo,BLMod=BLMod,
                                       method=optim.method,
                                       control = list(parscale = scale),
                                       hessian="T"))

        likelihood[i]<-mlest2$value*-1

      }

    }


    if(v==10){

      BLMod <- equation

      for(i in 1:length(sigh)){

        #sigh<-sigh[i]
        #likelihood<-vector()
        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ######################SUPPORT FUNCTIONS#########################

        nll_mef7<-function(pars,uplo,BLMod){


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
                         a=a,b=b,c=c,d=d,e=e,sigh=sigh[i],mux=mux,muy=muy,
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
        jdensup7<-function(X,BLMod,a,b,c,d,e,sigh=sigh[i],mux,muy,sdx,sdy,rcorr){
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
        coffcturb7<-function(x,mu,sig,A,C,sigh=sigh[i]){
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
        ###############################################################


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

        likelihood[i]<-mlest2$value*-1
      }
    }

  }

 ###############################

  if(plot==TRUE){
    plot(sigh,likelihood, xlab="Measurement error standard deviation",
       ylab="log-likelihood", main="Profile Likelihood", ...)}

 x<-list(likelihood=likelihood,Merror=sigh)
 x
}


