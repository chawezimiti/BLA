#' Predict boundary response
#'
#' This function predicts the most efficient response at a level of factor,
#' \code{x}, given the parameters of the fitted boundary line.
#'
#'
#' @param object An output in form of a list from the boundary line fitting using
#'   the \code{blqr()}, \code{blbin()}, \code{bolides()} or \code{cbvn()} functions.
#' @param x A numeric vector of values for the factor with which response is
#'   to be predicted.
#'
#' @returns A vector predicted value of response.
#'
#' @author Chawezi Miti <chawezi.miti@@nottingham.ac.uk>
#'
#' @export
#'
#' @examples
#'
#' x<-evapotranspiration$`ET(mm)`
#' y<-evapotranspiration$`yield(t/ha)`
#' z<-bolides(x,y, theta = c(0.5,0.02), model= "blm", xmax = 350)
#'
#' predictBL(z,x)
#'
#'
predictBL<-function(object,x){

  if(object$Model=="blm"){
    y<-tryCatch(lapply(x,
                       function(a,b) b$Parameters[1,1] + b$Parameters[2,1]*a,
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="lp"){
    y<-tryCatch(lapply(x,
                       function(a,b) min(b$Parameters[1,1] + b$Parameters[2,1]*a,b$Parameters[3,1],na.rm = F),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="logistic"|object[[1]]=="logistic"){

    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="inv-logistic"|object[[1]]=="inv-logistic"){

    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]- (b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x)))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="logisticfm"|object[[1]]=="logisticfm"){

    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]/(1+(b$Parameters[1,1]*exp(-b$Parameters[2,1]*x))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="double-logistic"|object[[1]]=="double-logistic"){

    y<-tryCatch(lapply(x,
            function(x,b) {
            (b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x))))-(b$Parameters[4,1]/(1+exp(b$Parameters[6,1]*(b$Parameters[5,1]-x))))
                         },
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="qd"){

    y<-tryCatch(lapply(x,
                       function(x,b) {b$Parameters[1,1] + b$Parameters[2,1]*x + b$Parameters[3,1]*x^2},
                       b=object),error=function(e) NA)
    return(unlist(y))
  }


  if(object$Model=="trapezium"){

    b<-object
    yr<-b$Parameters[1,1]+b$Parameters[2,1]*x
    yf<-b$Parameters[4,1]+b$Parameters[5,1]*x
    ym<-rep(b$Parameters[3,1],length(x))

    dat<-data.frame(yr,yf,ym)
    y<-apply(dat, 1, min)
    return(y)

  }

  if(object$Model=="schmidt"){

    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]-b$Parameters[1,1]*(x-b$Parameters[2,1])^2,
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="mit"){

    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[1,1]+b$Parameters[3,1]*(1-exp(-x/b$Parameters[2,1])),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }

  if(object$Model=="other"){

    y<-do.call(object$Equation, c(list(x=x),as.list(c(object$Parameters[,1]))))
    return(y)
  }

}
