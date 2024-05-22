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
#' @param start A numeric vector of initial starting values for optimization
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
#' @param sigh Standard deviation of the measurement error.
#' @param UpLo Selects the type of boundary. \code{"U"} fits the upper boundary and
#'   "L" fits the lower boundary.
#' @param model Selects the functional form of the boundary line. It includes
#'   \code{"blm"} for linear model, \code{"lp"} for linear plateau model, \code{"mit"}
#'   for the Mitscherlich model, \code{"schmidt"} for the Schmidt model,
#'   \code{"logistic"} for logistic model, \code{"logisticND"} for logistic model
#'   proposed by Nelder (1961), \code{"inv-logistic"} for the inverse logistic model,
#'   \code{"double-logistic"} for the double logistic model, \code{"qd"} for
#'   quadratic model and the \code{"trapezium"} for the trapezium model. For custom
#'   models, set \code{model = "other"}.
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
#'  where \eqn{\beta_1} is a scaling parameter, \eqn{\beta_2} is a
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
#' The common errors encountered due to poor start values \enumerate{
#' \item function cannot be evaluated at initial parameters
#' \item initial value in 'vmmin' is not finite}
#'
#'
#' @references
#'
#' Nelder, J.A. 1961. The fitting of a generalization of the logistic curve.
#' Biometrics 17: 89–110.
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
#' Phillips, B.F. & Campbell, N.A. 1968. A new method of fitting the von Bertelanffy
#' growth curve using data on the whelk. Dicathais, Growth 32: 317–329.
#'
#' Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line approach
#' to analyze N2O flux data from agricultural soils. Nutrient Cycling in Agroecosystems,
#' 57, 119-129.
#'
#' @author \enumerate{
#' \item Chawezi Miti \email{chawezi.miti@@nottingham.ac.uk}
#' \item Richard Murray Lark \email{murray.lark@@nottingham.ac.uk}
#' }
#' @import mvtnorm data.table numDeriv
#' @export
#'
#' @rdname cbvn
#' @usage
#' cbvn(vals,model="lp", equation=NULL, start, sigh, UpLo="U", optim.method="BFGS",
#'       Hessian=FALSE, plot=TRUE, line_smooth=1000, lwd=2, l_col="red",...)
#'
#' @examples
#'
#' x<-log(SoilP$P)
#' y<-SoilP$yield
#' vals<-data.frame(x,y)
#' start<-c(4,3,13.6, 35, -5,3,9,0.50,1.9,0.05)
#'
#' cbvn(vals,start=start,model = "trapezium", sigh = 0.7,
#'       xlab=expression("Phosphorus/ln(mg L"^-1*")"),
#'       ylab=expression("Yield/ t ha"^-1), pch=16,
#'       col="grey")
#'
cbvn<-function(vals, model="lp", equation=NULL, start, sigh, UpLo="U", optim.method="BFGS",
               Hessian=FALSE, plot=TRUE, line_smooth=1000, lwd=2, l_col="red",...){


  ########## Initial data preparations ##################################################

  vals <- na.omit(as.data.table(vals))

  sigh<-sigh # set value for the measurement error
  UpLo<-UpLo # set UpLo to "U" when fitting an upper boundary and "L" for a lower.

  ########## Fitting the three parameter model ###########################################

  if(model=="lp"|model=="mit"|model=="logistic"|model=="inv-logistic"|model=="schmidt"|model=="qd"|model=="logisticND"){

    if (length(start) != 8) stop("start must have exactly eight values")

    ## Define model functions-------------------------------------------------------------

      model_funcs <- list(
        lp = function(x, beta0, beta1, beta2) pmin(beta0, beta1 + beta2 * x),
        mit = function(x, beta0, beta1, beta2) beta0 - beta1 * beta2^x,
        logistic = function(x, beta0, beta1, beta2) beta0 / (1 + exp(beta2 * (beta1 - x))),
        `inv-logistic` = function(x, beta0, beta1, beta2) beta0 - (beta0 / (1 + exp(beta2 * (beta1 - x)))),
        logisticND = function(x, beta0, beta1, beta2) beta0 / (1 + beta1 * exp(-x * beta2)),
        schmidt = function(x, beta0, beta1, beta2) beta0 - beta1 * (x - beta2)^2,
        qd = function(x, beta0, beta1, beta2) beta1 + beta2 * x + beta0 * x^2
      )

      BLMod <- model_funcs[[model]]



    ## Define likelihood functions--------------------------------------------------------

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

    coffcturb <- function(x, mu, sig, a, c, sigh) {
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

    ## Optimization of the model----------------------------------------------------------

    mlest <- suppressWarnings(optim(start, nll_mef, uplo = UpLo, BLMod = BLMod,
                                    method = optim.method, hessian = TRUE))

    scale <- suppressWarnings(1 / abs(grad(nll_mef, mlest$par, uplo = UpLo, BLMod = BLMod)))

    mlest2 <- suppressWarnings(optim(mlest$par, nll_mef, uplo = UpLo, BLMod = BLMod,
                                     method = optim.method, control = list(parscale = scale),
                                     hessian = TRUE))

    AICbl <- (2 * 8) + (2 * mlest2$value)

    ## Extract parameters of the boundary line--------------------------------------------

    beta0 <- mlest2$par[3]
    beta1 <- mlest2$par[1]
    beta2 <- mlest2$par[2]

    ## Plotting the data for visualization------------------------------------------------

    if (plot) {
      plot(vals, ...)
      xdraw <- seq(min(vals$x), max(vals$x), length.out = line_smooth)
      ydraw <- BLMod(xdraw, beta0, beta1, beta2)
      lines(xdraw, ydraw, col = l_col, lwd = lwd)
    }

    ## Determine standard error for parameters--------------------------------------------

    estimates <- cbind(
      Estimate = mlest2$par,
      `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    )

    rownames(estimates) = c("\u03B2\u2081", "\u03B2\u2082", "\u03B2\u2080", "mux", "muy", "sdx", "sdy", "rcorr")


    ## Fitting the null model and compute its AIC-----------------------------------------

    start2 <- c(start[4:8], 0)
    mvnmlest <- suppressWarnings(optim(start2, nllmvn, method = optim.method))
    AImvn <- (2 * 5) + (2 * mvnmlest$value)

    ## Output preparation-----------------------------------------------------------------

    AikakeIC <- rbind(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn", "BL")
    colnames(AikakeIC)<-c("")

    equations <- list(
      lp = "y = min (\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080)",
      mit = "y = \u03B2\u2081 + \u03B2\u2080(1-exp(-x/\u03B2\u2081))",
      logistic = "y = \u03B2\u2080/1+[\u03B2\u2081exp(-\u03B2\u2082x)]",
      `inv-logistic` = "y = \u03B2\u2080/(1+exp(\u03B2\u2082(\u03B2\u2081-x)))",
      logisticND = "y = \u03B2\u2080/1+[\u03B2\u2081exp(-\u03B2\u2082*x)]",
      schmidt = "y = \u03B2\u2080 - \u03B2\u2081 (1-\u03B2\u2082)\u00B2)",
      qd = "y = \u03B2\u2081+\u03B2\u2082x+\u03B2\u2083x\u00B2"
    )
    Equation <- noquote(equations[[model]])

    result <- list(Model = model, Equation = Equation, Parameters = estimates, AIC = AikakeIC, Hessian = mlest2$hessian)
    class(result) <- "cm"
    return(result)

  }

  ########## Fitting the two parameter linear model ######################################

  if(model=="blm"){

    if (length(start) != 7) stop("start must have exactly seven values")

    ## Define model functions-------------------------------------------------------------


    blm <- function(x,beta0,beta1) beta0+beta1*x
    BLMod <- blm


    drawBL2<-function(x,beta0,beta1,BLMod){
      y<-sapply(x,BLMod,beta0=beta0,beta1=beta1)
      return(y)
    }


    ## Define likelihood functions--------------------------------------------------------

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

    coffcturb2 <- function(x, mu, sig, a, c, sigh) {
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


    ## Optimization of the model----------------------------------------------------------

    mlest <- suppressWarnings(optim(start, nll_mef2, uplo = UpLo, BLMod = BLMod,
                                    method = optim.method, hessian = TRUE))

    scale <- suppressWarnings(1 / abs(grad(nll_mef2, mlest$par, uplo = UpLo, BLMod = BLMod)))

    mlest2 <- suppressWarnings(optim(mlest$par, nll_mef2, uplo = UpLo, BLMod = BLMod,
                                     method = optim.method, control = list(parscale = scale),
                                     hessian = TRUE))

    ## Compute AIC value------------------------------------------------------------------

    AICbl<-(2*7)+(2*mlest2$value)

    ## Extract parameters of the boundary line-------------------------------------------

    beta0<-mlest2$par[1]
    beta1<-mlest2$par[2]

    ## Plotting the data for viewing------------------------------------------------------

    if(plot==TRUE){ plot(vals,...)

      x<-vals[,1]
      y<-vals[,2]

      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL2((xdraw),beta0,beta1,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)}

    ## Determine standard error for parameters--------------------------------------------

    estimates <- cbind(
      Estimate = mlest2$par,
      `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    )

    rownames(estimates) = c("\u03B2\u2081", "\u03B2\u2082", "mux", "muy", "sdx", "sdy", "rcorr")


    ##  Fitting the null model, an unbounded multivariate normal, and compute its AIC-----

    start2<-c(start[c(3,4,5,6)],0)
    mvnmlest<- suppressWarnings(optim( start2,nllmvn2, method=optim.method))
    AImvn<-(2*5)+(2*mvnmlest$value)


    ## Output Preparation-----------------------------------------------------------------

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")

    Equation<-noquote("y = \u03B2\u2081 + \u03B2\u2082x")

    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=mlest2$hessian)
    class(result)<-"cm"
    return(result)
  }


  ########## Fitting the five parameter trapezium model ##################################

  if(model=="trapezium"){

    if (length(start) != 10) stop("start must have exactly ten values")

    ## Define model functions------------------------------------------------------------

    trapezium <- function(x,beta0,beta1,beta2,beta3,beta4) pmin(beta0,beta1+beta2*x,beta3+beta4*x)
    BLMod <- trapezium

    drawBL3<-function(x,beta0,beta1,beta2,beta3,beta4,BLMod){
      y<-sapply(x,BLMod,beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
      return(y)
    }


    ## Define likelihood functions--------------------------------------------------------

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

    coffcturb3 <- function(x, mu, sig, a, c, sigh) {
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


    ## Optimization of the model---------------------------------------------------------

    mlest <- suppressWarnings(optim(start, nll_mef3, uplo = UpLo, BLMod = BLMod,
                                    method = optim.method, hessian = TRUE))

    scale <- suppressWarnings(1 / abs(grad(nll_mef3, mlest$par, uplo = UpLo, BLMod = BLMod)))

    mlest2 <- suppressWarnings(optim(mlest$par, nll_mef3, uplo = UpLo, BLMod = BLMod,
                                     method = optim.method, control = list(parscale = scale),
                                     hessian = TRUE))
    ## Determine the AIC-----------------------------------------------------------------

    AICbl<-(2*10)+(2*mlest2$value)

    ## Extract parameters of the boundary line-------------------------------------------

    beta0<-mlest2$par[3]
    beta1<-mlest2$par[1]
    beta2<-mlest2$par[2]
    beta3<-mlest2$par[4]
    beta4<-mlest2$par[5]

    ## Plotting the data for viewing-----------------------------------------------------

    if(plot==TRUE){ plot(vals,...)

      x<-vals[,1]
      y<-vals[,2]
      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL3((xdraw),beta0,beta1,beta2,beta3,beta4,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)
    }

    ## Determine standard error for parameters--------------------------------------------
    estimates <- cbind(
      Estimate = mlest2$par,
      `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    )

    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u2080","\u03B2\u2083","\u03B2\u2084","mux","muy","sdx","sdy","rcorr")

    ## Fitting the null model, an unbounded multivariate normal, and compute its AIC------

    start2<-c(start[c(6,7,8,9)],0)

    mvnmlest<-suppressWarnings(optim( start2,nllmvn3, method=optim.method))

    AImvn<-(2*5)+(2*mvnmlest$value)


    ## Output preparation---------------------------------------------------------------

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")

    Equation<-noquote("y = min (\u03B2\u2081 + \u03B2\u2082x, \u03B2\u2080, \u03B2\u2083 + \u03B2\u2084x)")
    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=mlest2$hessian)
    class(result)<-"cm"
    return(result)


  }

  ########## Fitting the six parameter model #############################################

  if(model=="double-logistic"){

    if (length(start) != 11) stop("start must have exactly eleven values")

    ## Define model functions-------------------------------------------------------------

    `double-logistic` <- function(x,beta01,beta02,beta1,beta2,beta3,beta4){
      (beta01/(1 + exp(beta2*(beta1-x)))) - (beta02/(1 + exp(beta4*(beta3-x))))
    }

    BLMod <- `double-logistic`

    drawBL4<-function(x,beta01,beta02,beta1,beta2,beta3,beta4,BLMod){
      y<-sapply(x,BLMod,beta01=beta01,beta02=beta02,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
      return(y)
    }

    ## Define likelihood functions-------------------------------------------------------

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

    coffcturb4 <- function(x, mu, sig, a, c, sigh) {
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


    ## Optimization of the model----------------------------------------------------------

    mlest <- suppressWarnings(optim(start, nll_mef4, uplo = UpLo, BLMod = BLMod,
                                    method = optim.method, hessian = TRUE))

    scale <- suppressWarnings(1 / abs(grad(nll_mef4, mlest$par, uplo = UpLo, BLMod = BLMod)))

    mlest2 <- suppressWarnings(optim(mlest$par, nll_mef4, uplo = UpLo, BLMod = BLMod,
                                     method = optim.method, control = list(parscale = scale),
                                     hessian = TRUE))

    ## Compute the Akaike Information Criterion for the fitted model----------------------

    AICbl<-(2*11)+(2*mlest2$value)

    ## Extract parameters of the boundary line--------------------------------------------

    beta1<-mlest2$par[1]
    beta2<-mlest2$par[2]
    beta01<-mlest2$par[3]
    beta02<-mlest2$par[4]
    beta3<-mlest2$par[5]
    beta4<-mlest2$par[6]

    ## Plotting the data for viewing------------------------------------------------------

    if(plot==TRUE){ plot(vals,...)
      x<-vals[,1]
      y<-vals[,2]
      xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
      ydraw<-drawBL4((xdraw),beta01,beta02,beta1,beta2,beta3,beta4,BLMod)
      lines(xdraw,ydraw,col=l_col,lwd=lwd)
    }

    # Determine standard error for parameters---------------------------------------------

    estimates <- cbind(
      Estimate = mlest2$par,
      `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)
    )
    rownames(estimates)<-c("\u03B2\u2081","\u03B2\u2082","\u03B2\u20801","\u03B2\u20802","\u03B2\u2083","\u03B2\u2084","mux","muy","sdx","sdy","rcorr")


    ## fitting the null model, an unbounded multivariate normal, and compute its AIC-----

    start2<-c(start[c(7,8,9,10)],0)
    mvnmlest<-suppressWarnings(optim( start2,nllmvn4, method=optim.method))
    AImvn<-(2*5)+(2*mvnmlest$value)


    ## Output preparation ----------------------------------------------------------------

    AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
    AikakeIC[,1]<-c(AImvn,AICbl)
    rownames(AikakeIC)<-c("mvn","BL")

    Equation<-noquote("y = {\u03B2\u20801/1+[exp(\u03B2\u2082*(\u03B2\u2081-x))]} - {\u03B2\u20801/1+[exp(\u03B2\u2084*(\u03B2\u2083-x))]} ")

    result<-list(Model=model,Equation=Equation, Parameters=estimates,AIC=AikakeIC, Hessian=mlest2$hessian)
    class(result)<-"cm"
    return(result)


  }

  ########## CUSTOM MODELS ###############################################################

  if(model=="other"){

    v<-length(start)
    if(length(start)>10) stop("Not set up for models with more than 5 parametrs. The argument start should contain less than 10 values")
    Equation<-equation # to print equation in output
    start<-unname(start) # removes names from start

    ### The three parameter models ######################################################

    if(v==8){

      BLMod <- equation

      drawBL5 <- function(x,a,b,c,BLMod){
        y<-sapply(x,BLMod,a=a,b=b,c=c)
        return(y)
      }

      ## Define likelihood functions------------------------------------------------------

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

      coffcturb5 <- function(x, mu, sig, A, C, sigh) {
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


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(start, nll_mef5, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef5, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef5, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      ## Compute the Akaike Information Criterion for the fitted model--------------------

      AICbl<-(2*8)+(2*mlest2$value)

      ## Extract parameters of the boundary line------------------------------------------

      a<-mlest2$par[1]
      b<-mlest2$par[2]
      c<-mlest2$par[3]


      ## Plotting the data for viewing----------------------------------------------------

      if(plot==TRUE){ plot(vals,...)
        x<-vals[,1]
        y<-vals[,2]
        xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
        ydraw<-drawBL5((xdraw),a,b,c,BLMod)
        lines(xdraw,ydraw,col=l_col,lwd=lwd)
      }

      ## Determine standard error for parameters------------------------------------------

      estimates <- cbind(
        Estimate = mlest2$par[1:3],
        `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[1:3]
      )
      rownames(estimates)<-c("a","b","c")

      distribution<-matrix(NA,5,2,dimnames=list(c(),c("Estimate","Standard error")))
      distribution[,1]<-mlest2$par[4:8]
      distribution[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[4:8]
      rownames(distribution)<-c("mux","muy","sdx","sdy","rcorr")

      ## Fitting the null model, an unbounded multivariate normal, and compute its AIC---

      start2<-c(start[c(4,5,6,7)],0)
      mvnmlest<-suppressWarnings(optim( start2,nllmvn5, method=optim.method))
      AImvn<-(2*5)+(2*mvnmlest$value) #

      ## Output preparation --------------------------------------------------------------

      AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
      AikakeIC[,1]<-c(AImvn,AICbl)
      rownames(AikakeIC)<-c("mvn","BL")

      result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Distribution=distribution, Hessian=mlest2$hessian)
      class(result)<-"cm"
      return(result)

    }


    #### Fitting the four parameter model ###############################################

    if(v==9){

      BLMod <- equation

      drawBL6 <- function(x,a,b,c,d,BLMod){
        y<-sapply(x,BLMod,a=a,b=b,c=c,d=d)
        return(y)
      }

      ## Define likelihood functions------------------------------------------------------

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

      coffcturb6 <- function(x, mu, sig, A, C, sigh) {
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


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(start, nll_mef6, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef6, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef6, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      ## Compute the Akaike Information Criterion for the fitted model--------------------

      AICbl<-(2*9)+(2*mlest2$value)

      ## Extract parameters of the boundary line-----------------------------------------

      a<-mlest2$par[1]
      b<-mlest2$par[2]
      c<-mlest2$par[3]
      d<-mlest2$par[4]


      ## Plotting the data for viewing----------------------------------------------------

      if(plot==TRUE){ plot(vals,...)
        x<-vals[,1]
        y<-vals[,2]
        xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
        ydraw<-drawBL6((xdraw),a,b,c,d,BLMod)
        lines(xdraw,ydraw,col=l_col,lwd=lwd)
      }

      ## Compute the standard error of the parameters-------------------------------------

      estimates <- cbind(
        Estimate = mlest2$par[1:4],
        `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[1:4]
      )
      rownames(estimates)<-c("a","b","c","d")

      distribution<-matrix(NA,5,2,dimnames=list(c(),c("Estimate","Standard error")))
      distribution[,1]<-mlest2$par[5:9]
      distribution[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[5:9]
      rownames(distribution)<-c("mux","muy","sdx","sdy","rcorr")

      ## Fitting the null model, an unbounded multivariate normal, and compute its AIC---

      start2<-c(start[c(4,5,6,7)],0)
      mvnmlest<-suppressWarnings(optim( start2,nllmvn6, method=optim.method))
      AImvn<-(2*5)+(2*mvnmlest$value)


      ## Output preparation --------------------------------------------------------------

      AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
      AikakeIC[,1]<-c(AImvn,AICbl)
      rownames(AikakeIC)<-c("mvn","BL")
      result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Distribution= distribution, Hessian=mlest2$hessian)
      class(result)<-"cm"
      return(result)

    }

    #### Fitting the five parameter model ################################################

    if(v==10){

      BLMod <- equation
      drawBL7<-function(x,a,b,c,d,e,BLMod){
        y<-sapply(x,BLMod,a=a,b=b,c=c,d=d,e=e)
        return(y)
      }


      ## Define likelihood functions------------------------------------------------------

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

      coffcturb7 <- function(x, mu, sig, A, C, sigh) {
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


      ## Optimization of the model--------------------------------------------------------

      mlest <- suppressWarnings(optim(start, nll_mef7, uplo = UpLo, BLMod = BLMod,
                                      method = optim.method, hessian = TRUE))

      scale <- suppressWarnings(1 / abs(grad(nll_mef7, mlest$par, uplo = UpLo, BLMod = BLMod)))

      mlest2 <- suppressWarnings(optim(mlest$par, nll_mef7, uplo = UpLo, BLMod = BLMod,
                                       method = optim.method, control = list(parscale = scale),
                                       hessian = TRUE))

      ## Compute the Akaike Information Criterion for the fitted model--------------------

      AICbl<-(2*10)+(2*mlest2$value)

      ## Extract parameters of the boundary line------------------------------------------

      a<-mlest2$par[1]
      b<-mlest2$par[2]
      c<-mlest2$par[3]
      d<-mlest2$par[4]
      e<-mlest2$par[5]

      ## Plotting the data for viewing----------------------------------------------------

      if(plot==TRUE){ plot(vals,...)

        x<-vals[,1]
        y<-vals[,2]
        xdraw=seq(min(x),max(x),(max(x)-min(x))/((max(x)-min(x))*line_smooth))
        ydraw<-drawBL7((xdraw),a,b,c,d,e,BLMod)
        lines(xdraw,ydraw,col=l_col,lwd=lwd)
      }

      ## Determine standard error for parameters------------------------------------------

      estimates <- cbind(
        Estimate = mlest2$par[1:5],
        `Standard error` = seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[1:5]
      )

      rownames(estimates)<-c("a","b","c","d","e")

      distribution<-matrix(NA,5,2,dimnames=list(c(),c("Estimate","Standard error")))
      distribution[,1]<-mlest2$par[6:10]
      distribution[,2]<-seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[6:10]
      rownames(distribution)<-c("mux","muy","sdx","sdy","rcorr")


      ## Fitting the null model, an unbounded multivariate normal, and compute its AIC----

      start2<-c(start[c(6,7,8,9)],0)

      mvnmlest<-suppressWarnings(optim( start2,nllmvn7,method=optim.method))

      AImvn<-(2*5)+(2*mvnmlest$value)


      ## Preparing the output-------------------------------------------------------------

      AikakeIC<-matrix(NA,2,1,dimnames=list(c(),c("")))
      AikakeIC[,1]<-c(AImvn,AICbl)
      rownames(AikakeIC)<-c("mvn","BL")

      result<-list(Model=model,Equation=equation, Parameters=estimates,AIC=AikakeIC, Distribution=distribution, Hessian=mlest2$hessian)
      class(result)<-"cm"
      return(result)


    }
  }
}
