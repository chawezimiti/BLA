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
#'   \item For the \code{"logistic"}, \code{"inv-logistic"} and \code{"logisticND"}
#'   models, it is a vector of length 8 arranged as scaling parameter, shape parameter,
#'   the maximum or plateau value, mean of \code{x}, mean of \code{y},
#'   standard deviation of \code{x}, standard deviation of \code{y} and the
#'   correlation of \code{x} and \code{y}.
#'   \item For the \code{"double-logistic"} model, it is a vector of length 11
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
#' @param sigh A vector of the suggested standard deviations of the measurement error
#'   values.
#' @param UpLo Selects the type of boundary. \code{"U"} fits the upper boundary and
#'   "L" fits the lower boundary.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for linear plateau model, \code{"mit"}
#'   for the Mitscherlich model, \code{"schmidt"} for the Schmidt model,
#'   \code{"logistic"} for logistic model, \code{"logisticND"} for logistic model
#'   proposed by Nelder (1961), \code{"inv-logistic"} for the inverse logistic
#'   model, \code{"double-logistic"} for the double logistic model, \code{"qd"} for
#'   quadratic model and the \code{"trapezium"} for the trapezium model. For custom
#'   models, set \code{model = "other"}.
#' @param equation A custom model function writen in the form of an R function. Applies
#'   only when argument \code{model="other"}, else it is \code{NULL}.
#' @param optim.method Describes the method used to optimize the model as in the
#'   \code{optim()} function. The methods include \code{"Nelder-Mead"}, \code{"BFGS"},
#'   \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} and \code{"Brent"}.
#' @param plot If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
#'   is not part of output (default is \code{TRUE}).
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
#' @import mvtnorm data.table numDeriv
#' @export
#' @rdname ble_profile
#' @usage
#' ble_profile(vals, sigh, model="lp", equation=NULL,  theta, UpLo="U",
#'              optim.method="BFGS", plot=TRUE)
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
ble_profile<-function(vals, sigh, model="lp", equation=NULL, theta, UpLo="U", optim.method="BFGS", plot=TRUE){

  cat("Note: This function may take longer to run.\n\n")

  ###### Initial data preparation ##################################

  vals <- na.omit(as.data.table(vals))

  UpLo=UpLo
  BLMod<-model
  likelihood<-vector()

  ### Fitting the three parameter model profile ##########################################

  if(model=="lp"|model=="mit"|model=="logistic"|model=="inv-logistic"|model=="logisticND"|model=="schmidt"|model=="qd"){

    if (length(theta) != 8) stop("theta must have exactly eight values")

    ## Define model functions------------------------------------------------------------------------------------

    model_funcs <- list(
      lp = function(x, beta0, beta1, beta2) pmin(beta0, beta1 + beta2 * x),
      mit = function(x, beta0, beta1, beta2) beta0 - beta1 * beta2^x,
      logistic = function(x, beta0, beta1, beta2) beta0 / (1 + exp(beta2 * (beta1 - x))),
      inv_logistic = function(x, beta0, beta1, beta2) beta0 - (beta0 / (1 + exp(beta2 * (beta1 - x)))),
      logisticND = function(x, beta0, beta1, beta2) beta0 / (1 + beta1 * exp(-x * beta2)),
      schmidt = function(x, beta0, beta1, beta2) beta0 - beta2 * (x - beta1)^2,
      qd = function(x, beta0, beta1, beta2) beta1 + beta2 * x + beta0 * x^2
    )

    BLMod <- model_funcs[[model]]


    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      ## Define likelihood functions------------------------------------------------------

      nll_mef <- function(pars, uplo, BLMod) {
        beta0 <- pars[3]
        beta1 <- pars[1]
        beta2 <- pars[2]
        mux <- pars[4]
        muy <- pars[5]
        sdx <- pars[6]
        sdy <- pars[7]
        rcorr <- pars[8]

        if (uplo == "U") {
          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy
          bet <- cov / (sdx^2)

          fx <- dnorm(vals[, x], mux, sdx)
          muyc <- muy + ((vals[, x] - mux) * bet)
          sdyc <- sdy * sqrt(1 - rho^2)

          c <- BLMod(vals[, x], beta0, beta1, beta2)
          fy_x <- coffcturb(vals[, y], muyc, sdyc, -Inf, c, sigh)

          fxy <- fy_x * fx
          -sum(log(fxy))
        } else {
          stop("Error, not set up for lower boundary")
        }
      }

      coffcturb <- function(x, mu, sig, a, c, sigh=sigh[i]) {
        k <- (mu - c) / sig
        d <- (mu - a) / sig
        alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
        beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
        gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))

        com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
        com2 <- gamma * exp(com1)
        f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))

        f <- f * pnorm(c, mu, sig)
        f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
      }

      nllmvn <- function(pars) {
        mux <- pars[1]
        muy <- pars[2]
        sdx <- pars[3]
        sdy <- pars[4]
        rcorr <- pars[5]

        rho <- tanh(rcorr)
        cov <- rho * sdx * sdy

        Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
        lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
        -sum(lliks)
      }

      ## Optimization of the model-------------------------------------------------------

      mlest <- suppressWarnings(optim(theta, nll_mef, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      likelihood[i]<-mlest2$value*-1

    }

  }

  ### Fitting the two parameter model profile ############################################

  if(model=="blm"){

    if (length(theta) != 7) stop("theta must have exactly seven values")

    ## Define model functions-------------------------------------------------------------

    blm <- function(x,beta0,beta1) beta0+beta1*x
    BLMod <-blm

    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      ## Define likelihood functions------------------------------------------------------

      nll_mef2 <- function(pars, uplo, BLMod) {
        beta0<-pars[1]
        beta1<-pars[2]
        mux<-pars[3]
        muy<-pars[4]
        sdx<-pars[5]
        sdy<-pars[6]
        rcorr<-pars[7]

        if (uplo == "U") {
          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy
          bet <- cov / (sdx^2)

          fx <- dnorm(vals[, x], mux, sdx)
          muyc <- muy + ((vals[, x] - mux) * bet)
          sdyc <- sdy * sqrt(1 - rho^2)

          c <- BLMod(vals[, x], beta0, beta1)
          fy_x <- coffcturb2(vals[, y], muyc, sdyc, -Inf, c, sigh)

          fxy <- fy_x * fx
          -sum(log(fxy))
        } else {
          stop("Error, not set up for lower boundary")
        }
      }

      coffcturb2 <- function(x, mu, sig, a, c, sigh=sigh[i]) {
        k <- (mu - c) / sig
        d <- (mu - a) / sig
        alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
        beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
        gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))

        com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
        com2 <- gamma * exp(com1)
        f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))

        f <- f * pnorm(c, mu, sig)
        f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
      }

      nllmvn2 <- function(pars) {
        mux <- pars[1]
        muy <- pars[2]
        sdx <- pars[3]
        sdy <- pars[4]
        rcorr <- pars[5]

        rho <- tanh(rcorr)
        cov <- rho * sdx * sdy

        Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
        lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
        -sum(lliks)
      }


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(theta, nll_mef2, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef2, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef2, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      likelihood[i]<-mlest2$value*-1

    }

  }

  ### Fitting the five parameter model ###################################################

  if(model=="trapezium"){

    if (length(theta) != 10) stop("theta must have exactly ten values")

    ## Define model functions-------------------------------------------------------------

    trapezium <- function(x,beta0,beta1,beta2,beta3,beta4) pmin(beta0,beta1+beta2*x,beta3+beta4*x)
    BLMod<-trapezium

    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals


      ## Define likelihood functions-----------------------------------------------------

      nll_mef3 <- function(pars, uplo, BLMod) {
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

        if (uplo == "U") {
          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy
          bet <- cov / (sdx^2)

          fx <- dnorm(vals[, x], mux, sdx)
          muyc <- muy + ((vals[, x] - mux) * bet)
          sdyc <- sdy * sqrt(1 - rho^2)

          c <- BLMod(vals[, x], beta0, beta1, beta2, beta3, beta4)
          fy_x <- coffcturb3(vals[, y], muyc, sdyc, -Inf, c, sigh)

          fxy <- fy_x * fx
          -sum(log(fxy))
        } else {
          stop("Error, not set up for lower boundary")
        }
      }

      coffcturb3 <- function(x, mu, sig, a, c, sigh=sigh[i]) {
        k <- (mu - c) / sig
        d <- (mu - a) / sig
        alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
        beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
        gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))

        com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
        com2 <- gamma * exp(com1)
        f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))

        f <- f * pnorm(c, mu, sig)
        f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
      }

      nllmvn3 <- function(pars) {
        mux <- pars[1]
        muy <- pars[2]
        sdx <- pars[3]
        sdy <- pars[4]
        rcorr <- pars[5]

        rho <- tanh(rcorr)
        cov <- rho * sdx * sdy

        Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
        lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
        -sum(lliks)
      }


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(theta, nll_mef3, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef3, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef3, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

    likelihood[i]<-mlest2$value*-1
    }

  }

  ### Fitting the six parameter double-logistic model ####################################

  if(model=="double-logistic"){


    if (length(theta) != 11) stop("theta must have exactly eleven values")

    ## Define model functions------------------------------------------------------------

    double_logistic <- function(x,beta01,beta02,beta1,beta2,beta3,beta4){
      (beta01/(1 + exp(beta2*(beta1-x)))) - (beta02/(1 + exp(beta4*(beta3-x))))
    }

    BLMod<-double_logistic


    for(i in 1:length(sigh)){

      theta<-theta
      UpLo=UpLo
      BLMod<-BLMod
      vals<-vals

      ## Define likelihood functions------------------------------------------------------

      nll_mef4 <- function(pars, uplo, BLMod) {
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

        if (uplo == "U") {
          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy
          bet <- cov / (sdx^2)

          fx <- dnorm(vals[, x], mux, sdx)
          muyc <- muy + ((vals[, x] - mux) * bet)
          sdyc <- sdy * sqrt(1 - rho^2)

          c <- BLMod(vals[, x], beta1, beta2, beta01, beta02, beta3, beta4)
          fy_x <- coffcturb4(vals[, y], muyc, sdyc, -Inf, c, sigh)

          fxy <- fy_x * fx
          -sum(log(fxy))
        } else {
          stop("Error, not set up for lower boundary")
        }
      }

      coffcturb4 <- function(x, mu, sig, a, c, sigh=sigh[i]) {
        k <- (mu - c) / sig
        d <- (mu - a) / sig
        alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
        beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
        gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))

        com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
        com2 <- gamma * exp(com1)
        f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))

        f <- f * pnorm(c, mu, sig)
        f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
      }

      nllmvn4 <- function(pars) {
        mux <- pars[1]
        muy <- pars[2]
        sdx <- pars[3]
        sdy <- pars[4]
        rcorr <- pars[5]

        rho <- tanh(rcorr)
        cov <- rho * sdx * sdy

        Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
        lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
        -sum(lliks)
      }


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(theta, nll_mef4, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef4, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef4, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      likelihood[i]<-mlest2$value*-1
    }

  }

  ### Fitting custom models ##############################################################

  if(model=="other"){

    v<-length(theta)
    if(v>10) stop("Not set up for models with more than 5 parametrs. The argument theta should contain less than 10 values")
    Equation<-equation # to print equation in output
    theta<-unname(theta) # removes names from theta

    if(v==8){

      ### The three parameter models -----------------------------------------------------

      BLMod <- equation

      for(i in 1:length(sigh)){

        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ## Define likelihood functions----------------------------------------------------

        nll_mef5 <- function(pars, uplo, BLMod) {
          a<-pars[1]
          b<-pars[2]
          c<-pars[3]
          mux<-pars[4]
          muy<-pars[5]
          sdx<-pars[6]
          sdy<-pars[7]
          rcorr<-pars[8]

          if (uplo == "U") {
            rho <- tanh(rcorr)
            cov <- rho * sdx * sdy
            bet <- cov / (sdx^2)

            fx <- dnorm(vals[, x], mux, sdx)
            muyc <- muy + ((vals[, x] - mux) * bet)
            sdyc <- sdy * sqrt(1 - rho^2)

            C <- BLMod(vals[, x], a,b,c)
            fy_x <- coffcturb5(vals[, y], muyc, sdyc, -Inf, C, sigh)

            fxy <- fy_x * fx
            -sum(log(fxy))
          } else {
            stop("Error, not set up for lower boundary")
          }
        }

        coffcturb5 <- function(x, mu, sig, A, C, sigh=sigh[i]) {
          k <- (mu - C) / sig
          D <- (mu - A) / sig
          alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
          beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
          gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(D) - pnorm(k)))

          com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
          com2 <- gamma * exp(com1)
          f <- com2 * (pnorm((x - A - alpha) / beta) - pnorm((x - C - alpha) / beta))

          f <- f * pnorm(C, mu, sig)
          f + (dnorm((x - C), 0, sigh) * (1 - pnorm(C, mu, sig)))
        }

        nllmvn5 <- function(pars) {
          mux <- pars[1]
          muy <- pars[2]
          sdx <- pars[3]
          sdy <- pars[4]
          rcorr <- pars[5]

          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy

          Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
          lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
          -sum(lliks)
        }


        ## Optimization of the model------------------------------------------------------

        mlest <- suppressWarnings(optim(theta, nll_mef5, uplo = UpLo, BLMod = BLMod,
                                        method = optim.method, hessian = TRUE))

        scale <- suppressWarnings(1 / abs(grad(nll_mef5, mlest$par, uplo = UpLo, BLMod = BLMod)))

        mlest2 <- suppressWarnings(optim(mlest$par, nll_mef5, uplo = UpLo, BLMod = BLMod,
                                         method = optim.method, control = list(parscale = scale),
                                         hessian = TRUE))

        likelihood[i]<-mlest2$value*-1

      }
    }

    ### Fitting the four parameter model ------------------------------------------------

    if(v==9){

      BLMod <- equation

      for(i in 1:length(sigh)){
        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ## Define likelihood functions----------------------------------------------------

        nll_mef6 <- function(pars, uplo, BLMod) {
          a<-pars[1]
          b<-pars[2]
          c<-pars[3]
          d<-pars[4]
          mux<-pars[5]
          muy<-pars[6]
          sdx<-pars[7]
          sdy<-pars[8]
          rcorr<-pars[9]

          if (uplo == "U") {
            rho <- tanh(rcorr)
            cov <- rho * sdx * sdy
            bet <- cov / (sdx^2)

            fx <- dnorm(vals[, x], mux, sdx)
            muyc <- muy + ((vals[, x] - mux) * bet)
            sdyc <- sdy * sqrt(1 - rho^2)

            C <- BLMod(vals[, x], a,b,c,d)
            fy_x <- coffcturb6(vals[, y], muyc, sdyc, -Inf, C, sigh)

            fxy <- fy_x * fx
            -sum(log(fxy))
          } else {
            stop("Error, not set up for lower boundary")
          }
        }

        coffcturb6 <- function(x, mu, sig, A, C, sigh=sigh[i]) {
          k <- (mu - C) / sig
          D <- (mu - A) / sig
          alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
          beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
          gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(D) - pnorm(k)))

          com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
          com2 <- gamma * exp(com1)
          f <- com2 * (pnorm((x - A - alpha) / beta) - pnorm((x - C - alpha) / beta))

          f <- f * pnorm(C, mu, sig)
          f + (dnorm((x - C), 0, sigh) * (1 - pnorm(C, mu, sig)))
        }

        nllmvn6 <- function(pars) {
          mux <- pars[1]
          muy <- pars[2]
          sdx <- pars[3]
          sdy <- pars[4]
          rcorr <- pars[5]

          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy

          Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
          lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
          -sum(lliks)
        }


        ## Optimization of the model------------------------------------------------------

        mlest <- suppressWarnings(optim(theta, nll_mef6, uplo = UpLo, BLMod = BLMod,
                                        method = optim.method, hessian = TRUE))

        scale <- suppressWarnings(1 / abs(grad(nll_mef6, mlest$par, uplo = UpLo, BLMod = BLMod)))

        mlest2 <- suppressWarnings(optim(mlest$par, nll_mef6, uplo = UpLo, BLMod = BLMod,
                                         method = optim.method, control = list(parscale = scale),
                                         hessian = TRUE))

        likelihood[i]<-mlest2$value*-1

      }

    }

    ### Fitting the five parameter model ------------------------------------------------

    if(v==10){

      BLMod <- equation

      for(i in 1:length(sigh)){

        theta<-theta
        UpLo=UpLo
        BLMod<-BLMod
        vals<-vals

        ## Define likelihood functions----------------------------------------------------

        nll_mef7 <- function(pars, uplo, BLMod) {
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

          if (uplo == "U") {
            rho <- tanh(rcorr)
            cov <- rho * sdx * sdy
            bet <- cov / (sdx^2)

            fx <- dnorm(vals[, x], mux, sdx)
            muyc <- muy + ((vals[, x] - mux) * bet)
            sdyc <- sdy * sqrt(1 - rho^2)

            C <- BLMod(vals[, x], a,b,c,d,e)
            fy_x <- coffcturb7(vals[, y], muyc, sdyc, -Inf, C, sigh)

            fxy <- fy_x * fx
            -sum(log(fxy))
          } else {
            stop("Error, not set up for lower boundary")
          }
        }

        coffcturb7 <- function(x, mu, sig, A, C, sigh=sigh[i]) {
          k <- (mu - C) / sig
          D <- (mu - A) / sig
          alpha <- (sigh^2 * (x - mu)) / (sigh^2 + sig^2)
          beta <- sqrt((sigh^2 * sig^2) / (sigh^2 + sig^2))
          gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(D) - pnorm(k)))

          com1 <- -((x - mu)^2) / (2 * (sigh^2 + sig^2))
          com2 <- gamma * exp(com1)
          f <- com2 * (pnorm((x - A - alpha) / beta) - pnorm((x - C - alpha) / beta))

          f <- f * pnorm(C, mu, sig)
          f + (dnorm((x - C), 0, sigh) * (1 - pnorm(C, mu, sig)))
        }

        nllmvn7 <- function(pars) {
          mux <- pars[1]
          muy <- pars[2]
          sdx <- pars[3]
          sdy <- pars[4]
          rcorr <- pars[5]

          rho <- tanh(rcorr)
          cov <- rho * sdx * sdy

          Sigma <- matrix(c(sdx^2, cov, cov, sdy^2), 2, 2)
          lliks <- dmvnorm(vals, mean = c(mux, muy), sigma = Sigma, log = TRUE)
          -sum(lliks)
        }


        ## Optimization of the model------------------------------------------------------

        mlest <- suppressWarnings(optim(theta, nll_mef7, uplo = UpLo, BLMod = BLMod,
                                        method = optim.method, hessian = TRUE))

        scale <- suppressWarnings(1 / abs(grad(nll_mef7, mlest$par, uplo = UpLo, BLMod = BLMod)))

        mlest2 <- suppressWarnings(optim(mlest$par, nll_mef7, uplo = UpLo, BLMod = BLMod,
                                         method = optim.method, control = list(parscale = scale),
                                         hessian = TRUE))

        likelihood[i]<-mlest2$value*-1
      }
    }

  }

 ### Output preparation ##################################################################

  if(plot==TRUE){
    plot(sigh,likelihood, xlab="Measurement error standard deviation",
       ylab="log-likelihood", main="Profile Likelihood")}

 x<-list(likelihood=likelihood,Merror=sigh)
 x
}


