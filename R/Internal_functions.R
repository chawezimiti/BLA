#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@param uplo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'
#'@returns Parameters of nll_mef
#' @keywords internal
#'
#'@author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
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
beta0<-pars[1]
beta1<-pars[2]
beta2<-pars[3]
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x input in nll_mef function
#'@param UpLo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'@returns Parameters of par_nll_mef
#' @keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param X input in function
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'
#'@returns Parameters of jdensup
#'@keywords internal
jdensup<-function(X,BLMod,beta0,beta1,beta2,sigh,mux,muy,sdx,sdy,rcorr){
#########################################################################

# joint density of observed values x and y given beta0,beta1,beta2 as bl
# parameters (plateau, intercept and slope of bounded linear model).
# sigh ss measurement error
# mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
# NB parameterization of rcorr to keep in [-1,1]

x<-X[1]
y<-X[2]

BLGen<-match.fun(BLMod)

rho<-tanh(rcorr)
cov<-rho*sdx*sdy
bet<-cov/(sdx*sdx)

fx<-dnorm(x,mux,sdx)

muyc<-muy+((x-mux)*bet)
sdyc<-sdy*sqrt(1-(rho*rho))

c<-BLGen(x,beta0,beta1,beta2)

fy_x<-coffcturb(y,muyc,sdyc,-Inf,c,sigh)

fxy<-fy_x*fx
return(log(fxy))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@export
#'@returns Parameters of the model
#'@keywords internal
lp<-function(x,beta0,beta1,beta2){
#########################################################################
return(min(beta0,beta1+beta2*x))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@export
#'@returns Parameters of the model.
#'@keywords internal
mit<-function(x,beta0,beta1,beta2){
#########################################################################
return(beta1+beta0*(1-exp(-x/beta2)))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@export
#'@returns Parameters of the model
#'@keywords internal
blm<-function(x,beta0,beta1){
  #########################################################################
  return(beta0+beta1*x)
}

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 is a parameter describing the possible maximum response
#'@param beta1 is value of x at which the response is maximum
#'@param beta2 is a scaling parameter, which is zero if the boundary is a constant
#'@export
#'@returns Parameters of the model.
#'@keywords internal
schmidt<-function(x,beta0,beta1,beta2){
  #########################################################################
  return(beta0-beta1*(x-beta2)*(x-beta2))
}

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 is a parameter describing the possible maximum response
#'@param beta1 is value of x at which the response is maximum
#'@param beta2 is a scaling parameter, which is zero if the boundary is a constant
#'@export
#'@returns Parameters of the model.
#'@keywords internal
qd<-function(x,beta0,beta1,beta2){
  #########################################################################
  return(beta1+beta2*x+beta0*x*x)
}

######################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 is a parameter describing the possible maximum response
#'@param beta1 is a scaling parameter
#'@param beta2 is a parameter describing the shape of the model
#'@export
#'@returns Parameters of the model.
#'@keywords internal

BL_logistic<-function(x,beta0,beta1,beta2){
  #########################################################################
  return(beta0/(1+exp(beta2*(beta1-x))))
}

######################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 is a parameter describing the possible maximum response
#'@param beta1 is a scaling parameter
#'@param beta2 is a parameter describing the shape of the model
#'@export
#'@returns Parameters of the model.
#'@keywords internal

BL_inv_logistic<-function(x,beta0,beta1,beta2){
  #########################################################################
  return(beta0-(beta0/(1+exp(beta2*(beta1-x)))))
}
#########################################################################

######################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 is a parameter describing the possible maximum response
#'@param beta1 is a scaling parameter
#'@param beta2 is a parameter describing the shape of the model
#'@export
#'@returns Parameters of the model.
#'@keywords internal

BL_logisticfm<-function(x,beta0,beta1,beta2){
  #########################################################################
  return(beta0/(1+beta1*exp(-x*beta2)))
}
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'
#'@returns draws the model on plot
#'@keywords internal
drawBL<-function(x,beta0,beta1,beta2,BLMod){
#########################################################################
BLGen<-match.fun(BLMod)
y<-sapply(x,BLGen,beta0=beta0,beta1=beta1,beta2=beta2)
return(y)
}
########################################################################
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param sigh measurement error
#'@param mu mean
#'@param sig fitting parameter
#'@param a fitting parameter
#'@param c fitting parameter
#'
#'@returns Parameters of coffcturb
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'
#'@returns Parameters of null model
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'
#'@returns Parameters of the max response model
#'@keywords internal
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

#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param X input in function
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'@param ymax maximumresponse
#'
#'@returns Parameters of jdens_maxyield
#'@keywords internal
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
#########################################################################

######Additions from linear model######################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@param uplo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'
#'@returns Parameters nll_mef2 for linear model
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param UpLo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'@param x independent variable
#'
#'@returns Parameters of par_nll_mef2
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param X input in function
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'
#'@returns Parameters of jdensup2
#'@keywords internal
jdensup2<-function(X,BLMod,beta0,beta1,sigh,mux,muy,sdx,sdy,rcorr){
  #########################################################################

  # joint density of observed values x and y given beta0,beta1,beta2 as bl
  # parameters (plateau, intercept and slope of bounded linear model).
  # sigh ss measurement error
  # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
  # NB parameterization of rcorr to keep in [-1,1]

  x<-X[1]
  y<-X[2]

  BLGen2<-match.fun(BLMod)

  rho<-tanh(rcorr)
  cov<-rho*sdx*sdy
  bet<-cov/(sdx*sdx)

  fx<-dnorm(x,mux,sdx)

  muyc<-muy+((x-mux)*bet)
  sdyc<-sdy*sqrt(1-(rho*rho))

  c<-BLGen2(x,beta0,beta1)

  fy_x<-coffcturb(y,muyc,sdyc,-Inf,c,sigh)

  fxy<-fy_x*fx
  return(log(fxy))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'
#'@returns Parameters of BL_lm.
#'@keywords internal
blm<-function(x,beta0,beta1){
  #########################################################################
  return(beta0+beta1*x)
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'
#'@returns draws model line on plot
#'@keywords internal
drawBL2<-function(x,beta0,beta1,BLMod){
  #########################################################################
  BLGen2<-match.fun(BLMod)
  y<-sapply(x,BLGen2,beta0=beta0,beta1=beta1)
  return(y)
}


############## Trapezium ###################################################

#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@param uplo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'
#'@returns Parameters of nll_mef
#'@keywords internal
#'@author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
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
  beta0<-pars[1]
  beta1<-pars[2]
  beta2<-pars[3]
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x input in nll_mef3 function
#'@param UpLo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'@returns Parameters of par_nll_mef
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param X input in function
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'
#'@returns Parameters of jdensup
#'@keywords internal
jdensup3<-function(X,BLMod,beta0,beta1,beta2,beta3,beta4,sigh,mux,muy,sdx,sdy,rcorr){
  #########################################################################

  # joint density of observed values x and y given beta0,beta1,beta2 as bl
  # parameters (plateau, intercept and slope of bounded linear model).
  # sigh ss measurement error
  # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
  # NB parameterization of rcorr to keep in [-1,1]

  x<-X[1]
  y<-X[2]

  BLGen3<-match.fun(BLMod)

  rho<-tanh(rcorr)
  cov<-rho*sdx*sdy
  bet<-cov/(sdx*sdx)

  fx<-dnorm(x,mux,sdx)

  muyc<-muy+((x-mux)*bet)
  sdyc<-sdy*sqrt(1-(rho*rho))

  c<-BLGen3(x,beta0,beta1,beta2,beta3,beta4)

  fy_x<-coffcturb3(y,muyc,sdyc,-Inf,c,sigh)

  fxy<-fy_x*fx
  return(log(fxy))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@keywords internal
#'@export
#'@returns Parameters of the model
trapezium<-function(x,beta0,beta1,beta2,beta3,beta4){
  #########################################################################
  return(min(beta0,beta1+beta2*x,beta3+beta4*x ))
}
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param x independent variable
#'@param beta0 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@keywords internal
#'
#'@returns draws the model on plot
drawBL3<-function(x,beta0,beta1,beta2,beta3,beta4,BLMod){
  #########################################################################
  BLGen3<-match.fun(BLMod)
  y<-sapply(x,BLGen3,beta0=beta0,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
  return(y)
}
########################################################################
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param sigh measurement error
#'@param mu mean
#'@param sig fitting parameter
#'@param a fitting parameter
#'@param c fitting parameter
#'@keywords internal
#'
#'@returns Parameters of coffcturb
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@keywords internal
#'
#'@returns Parameters of null model
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
  nllmvn3<--sum(lliks)

  return(nllmvn3)
}

########################################################################


#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@keywords internal
#'
#'@returns Parameters of the max response model
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

#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param X input in function
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'@param ymax maximumresponse
#'@keywords internal
#'
#'@returns Parameters of jdens_maxyield
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
#########################################################################

############## Double-logistic###################################################

#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@param uplo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'
#'@returns Parameters of nll_mef
#'@keywords internal
#'@author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x input in nll_mef3 function
#'@param UpLo determines whether it is upper or lower boundary
#'@param BLMod determines the boundary line model fitted
#'@returns Parameters of par_nll_mef
#'@keywords internal
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param X input in function
#'@param beta01 model parameter
#'@param beta02 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'
#'@returns Parameters of jdensup
#'@keywords internal
jdensup4<-function(X,BLMod,beta01,beta02,beta1,beta2,beta3,beta4,sigh,mux,muy,sdx,sdy,rcorr){
  #########################################################################

  # joint density of observed values x and y given beta0,beta1,beta2 as bl
  # parameters (plateau, intercept and slope of bounded linear model).
  # sigh ss measurement error
  # mux,muy,sdx,sdy,rcorr as parameters of underlying bivariate normal rv
  # NB parameterization of rcorr to keep in [-1,1]

  x<-X[1]
  y<-X[2]

  BLGen4<-match.fun(BLMod)

  rho<-tanh(rcorr)
  cov<-rho*sdx*sdy
  bet<-cov/(sdx*sdx)

  fx<-dnorm(x,mux,sdx)

  muyc<-muy+((x-mux)*bet)
  sdyc<-sdy*sqrt(1-(rho*rho))

  c<-BLGen4(x,beta01,beta02,beta1,beta2,beta3,beta4)

  fy_x<-coffcturb4(y,muyc,sdyc,-Inf,c,sigh)

  fxy<-fy_x*fx
  return(log(fxy))
}
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param beta01 model parameter
#'@param beta02 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@keywords internal
#'@export
#'@returns Parameters of the model
BL_double_logistic<-function(x,beta01,beta02,beta1,beta2,beta3,beta4){
  #########################################################################

  return((beta01/(1 + exp(beta2*(beta1-x)))) - (beta02/(1 + exp(beta4*(beta3-x)))))
}
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param BLMod determines the boundary line model fitted
#'@param x independent variable
#'@param beta01 model parameter
#'@param beta02 model parameter
#'@param beta1 model parameter
#'@param beta2 model parameter
#'@param beta3 model parameter
#'@param beta4 model parameter
#'@keywords internal
#'
#'@returns draws the model on plot
drawBL4<-function(x,beta01,beta02,beta1,beta2,beta3,beta4,BLMod){
  #########################################################################
  BLGen4<-match.fun(BLMod)
  y<-sapply(x,BLGen4,beta01=beta01,beta02=beta02,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4)
  return(y)
}
########################################################################
#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param x independent variable
#'@param sigh measurement error
#'@param mu mean
#'@param sig fitting parameter
#'@param a fitting parameter
#'@param c fitting parameter
#'@keywords internal
#'
#'@returns Parameters of coffcturb
coffcturb4<-function(x,mu,sig,a,c,sigh){
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
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@keywords internal
#'
#'@returns Parameters of null model
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


#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param pars input in nll_mef function
#'@keywords internal
#'
#'@returns Parameters of the max response model
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

#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param X input in function
#'@param sigh measurement error
#'@param mux mean of x
#'@param muy mean of y
#'@param sdx standard deviation of x
#'@param sdy standard deviation of y
#'@param rcorr correlation of x and y
#'@param ymax maximumresponse
#'@keywords internal
#'
#'@returns Parameters of jdens_maxyield
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
#########################################################################

#########################################################################
#' Internal functions
#'
#'This is a set of functions that support the censored bivariate normal model.
#'
#'@param hessian If true hessian is used
#'@param silent condtion of matrix
#'@param a determines the boundary line model fitted
#'@keywords internal
#'
seHessian<-function(a, hessian = FALSE, silent = FALSE){
  namesp <- colnames(a)
  mathessian <- a
  mathessian <- ifelse(mathessian == -Inf, -1e+09, mathessian)
  mathessian <- ifelse(mathessian == +Inf, 1e+09, mathessian)
  sigma <- try(solve(mathessian), silent = TRUE)
  if (inherits(sigma, "try-error")) {
    if (!silent)
      warning("Error in Hessian matrix inversion")
    mathessianx <- try(as.matrix(getFromNamespace("nearPD",
                                                  ns = "Matrix")(mathessian)$mat), silent = TRUE)
    if (inherits(mathessianx, "try-error")) {
      if (!silent)
        warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
      sigma <- try(ginv(mathessian), silent = TRUE)
      if (is.null(colnames(sigma)) | is.null(rownames(sigma))) {
        colnames(sigma) <- rownames(sigma) <- colnames(mathessian)
      }
    }
    else {
      if (!silent)
        warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
      sigma <- try(solve(mathessianx), silent = TRUE)
    }
  }
  if (!inherits(sigma, "try-error")) {
    if (all(diag(sigma) >= 0)) {
      res <- sqrt(diag(sigma))
    }
    else {
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      res <- structure((matrix(rep(1, nrow(R)), nrow = 1,
                               byrow = TRUE) %*% R)[1, ], .Names = colnames(mathessian))
      if (any(res < 0)) {
        d <- diag(as.matrix(getFromNamespace("nearPD",
                                             ns = "Matrix")(sigma)$mat))
        names(d) <- colnames(mathessian)
        res <- ifelse(d < 0, NA, sqrt(d))
      }
      if (any(is.na(res))) {
        a <- sigma
        n = dim(a)[1]
        root = matrix(0, n, n)
        for (i in 1:n) {
          sum = 0
          if (i > 1) {
            sum = sum(root[i, 1:(i - 1)]^2)
          }
          x = a[i, i] - sum
          if (x < 0) {
            x = 0
          }
          root[i, i] = sqrt(x)
          if (i < n) {
            for (j in (i + 1):n) {
              if (root[i, i] == 0) {
                x = 0
              }
              else {
                sum = 0
                if (i > 1) {
                  sum = root[i, 1:(i - 1)] %*% t(t(root[j,
                                                        1:(i - 1)]))
                }
                x = (a[i, j] - sum)/root[i, i]
              }
              root[j, i] = x
            }
          }
        }
        colnames(root) <- rownames(root) <- colnames(mathessian)
        pseudoV <- root %*% t(root)
        d <- diag(pseudoV)
        if (any(d != 0) & all(d >= 0)) {
          res <- sqrt(d)
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
        }
        else {
          if (!silent)
            warning("Approximation of Cholesky matrix based on Rebonato and Jackel (2000)")
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
          newMat <- a
          cholError <- TRUE
          iter <- 0
          while (cholError) {
            iter <- iter + 1
            newEig <- eigen(newMat)
            newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
            newMat <- newEig$vectors %*% diag(newEig2) %*%
              t(newEig$vectors)
            newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
            cholStatus <- try(u <- chol(newMat), silent = TRUE)
            cholError <- ifelse(inherits(cholStatus,
                                         "try-error"), TRUE, FALSE)
          }
          root <- cholStatus
          colnames(root) <- rownames(root) <- colnames(mathessian)
          pseudoV <- root %*% t(root)
          res <- sqrt(diag(pseudoV))
        }
      }
    }
  }
  SEInf <- namesp[!namesp %in% names(res)]
  res <- c(res, structure(rep(+Inf, length(SEInf)), .Names = SEInf))
  if (hessian) {
    return(list(SE = res, hessian = mathessian))
  }
  else {
    return(res)
  }
}


