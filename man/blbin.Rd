% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Binning.R
\name{blbin}
\alias{blbin}
\title{Binning method for determining the boundary line model}
\usage{
blbin(x,y,bins, model="explore", equation=NULL, start, tau=0.95,
      optim.method="Nelder-Mead", xmin=min(bound$x), xmax=max(bound$x),plot=TRUE,
      bp_col="red", bp_pch=16, bl_col="red", lwd=1,line_smooth=1000,...)
}
\arguments{
\item{x}{A numeric vector of values for the independent variable.}

\item{y}{A numeric vector of values for the response variable.}

\item{bins}{A numeric vector of length 3 or 4 that determines the size of sections.
The first and second values give the range of the data to be binned while the
third and fourth values give the width of the bins and the step size
respectively. If only three values are provided, the step size is assumed to be
equal to bin width.}

\item{model}{Selects the functional form of the boundary line. It includes
\code{"explore"} as default, \code{"blm"} for linear model, \code{"lp"} for
linear plateau model, \code{"mit"} for the Mitscherlich model, \code{"schmidt"}
for the Schmidt model, \code{"logistic"} for logistic model, \code{"logisticND"}
for logistic model proposed by Nelder (1961), \code{"inv-logistic"} for
the inverse logistic model, \code{"double-logistic"} for the double logistic model,
\code{"qd"} for quadratic model and the \code{"trapezium"} for the trapezium model.
The \code{"explore"} is used to check the position of boundary points in each bin
so that the correct \code{model} can be applied. For custom models, set
\code{model = "other"}.}

\item{equation}{A custom model function writen in the form of an R function. Applies
only when argument \code{model="other"}, else it is \code{NULL}.}

\item{start}{A numeric vector of initial starting values for optimization
in fitting the boundary model. Its length and arrangement depend on the
suggested model: \itemize{
\item For the \code{"blm"} model, it is a vector of length 2 arranged as intercept
and slope.
\item For the \code{"lp"} model, it is a vector of length 3 arranged as intercept,
slope and maximum response.
\item For the \code{"logistic"} and \code{"inv-logistic"} models, it is a
vector of length 3 arranged as the scaling parameter, shape parameter and maximum
response.
\item For the \code{"logisticND"} model proposed by Nelder (1961), it is a
vector of length 3 arranged as the scaling parameter, shape parameter and maximum
response.
\item For the \code{"double-logistic"} model, it is a vector of length 6 arranged
as the scaling parameter one, shape parameter one, maximum response, maximum
response, scaling parameter two and shape parameter two.
\item For the \code{"qd"} model, it is a vector of length 3 arranged as constant,
linear coefficient and quadratic coefficient.
\item For the \code{"trapezium"} model, it is a vector of length 3  arranged as
intercept one, slope one, maximum response, intercept two and slope two.
\item For the \code{"mit"} model, it is a vector of length 3 arranged as the
intercept, shape parameter and the maximum response.
\item For the \code{"schmidt"} model, it is a vector of length 3 arranged as scaling
parameter, shape parameter (x-value at maximum response ) and maximum response.}}

\item{tau}{A percentile value (0-1) that represents the boundary point within each
bin (default is \code{tau = 0.95}).}

\item{optim.method}{Describes the method used to optimize the model as in the
\code{optim()} function. The methods include \code{"Nelder-Mead"}, \code{"BFGS"},
\code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} and \code{"Brent"}.}

\item{xmin}{Numeric value that describes the minimum \code{x} value to which the
boundary line is to be fitted (default is \code{min(x)}).}

\item{xmax}{A numeric value that describes the maximum \code{x} value to which the
boundary line is to be fitted (default is \code{max(x)}). \code{xmin} and
\code{xmax} determine the subset of the data set used to fit boundary model.}

\item{plot}{If \code{TRUE}, a plot is part of the output. If \code{FALSE}, plot
is not part of output (default is \code{TRUE}).}

\item{bp_col}{Selects the color of the boundary points.}

\item{bp_pch}{Point character as \code{pch} of the \code{plot()} function. It controls
the shape of the boundary points on plot (\code{bp_pch = 16} as default).}

\item{bl_col}{Colour of the boundary line.}

\item{lwd}{Determines the thickness of the boundary line on the plot (default is 1).}

\item{line_smooth}{Parameter that describes the smoothness of the boundary line.
(default is 1000). The higher the value, the smoother the line.}

\item{...}{Additional graphical parameters as in the \code{par()} function.}
}
\value{
A list of length 5 consisting of the fitted model, equation form, parameters
  of the boundary line, the residue mean square and the boundary points. Additionally,
  a graphical representation of the boundary line on the scatter plot is produced.
}
\description{
This function fits a boundary model to the upper bounds of a scatter plot of
\code{x} and \code{y} based on the binning method. The data are first divided
into equal sized sections in the x-axis and a boundary point in each section
is selected based on a set criteria (e.g. 0.90, 0.95 or 0.99 percentile of
\code{y} among other criteria). A model is then fitted to the resulting boundary
points by the least squares method. This is done using optimization procedure
and hence requires some starting guess parameters for the proposed model.
}
\details{
Some inbuilt models are available for the \code{blbin()} function. The
\code{"explore"} option for the argument \code{model} generates a plot showing the
location of the boundary points selected by the binning procedure. This helps to
identify which model type is suitable to fit as a boundary line. The suggest model
forms are as follows: \enumerate{
 \item Linear model (\code{"blm"})
 \deqn{y=\beta_1 + \beta_2x}
 where \eqn{\beta_1} is the intercept and \eqn{\beta_2} is the slope.

 \item Linear plateau model (\code{"lp"})
 \deqn{y= {\rm min}(\beta_1+\beta_2x, \beta_0)}
 where \eqn{\beta_1} is the intercept , \eqn{\beta_2} is the slope  and \eqn{\beta_0}
 is the maximum response.

 \item The logistic (\code{"logistic"}) and inverse logistic (\code{"inv-logistic"})
 models
 \deqn{ y= \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}}
 \deqn{ y= \beta_0 - \frac{\beta_0}{1+e^{\beta_2(\beta_1-x)}}}
 where \eqn{\beta_1} is a scaling parameter , \eqn{\beta_2} is a shape parameter
 and \eqn{\beta_0} is the maximum response.

 \item Logistic model (\code{"logisticND"})  (Nelder (2009))
 \deqn{ y= \frac{\beta_0}{1+(\beta_1 \times e^{-\beta_2x})}}
  where \eqn{\beta_1} is a scaling parameter, \eqn{\beta_2} is a shape
  parameter and \eqn{\beta_0} is the maximum response.

 \item Double logistic model (\code{"double-logistic"})
 \deqn{ y= \frac{\beta_{0,1}}{1+e^{\beta_2(\beta_1-x)}} -
 \frac{\beta_{0,2}}{1+e^{\beta_4(\beta_3-x)}}}
 where \eqn{\beta_1} is a scaling parameter one, \eqn{\beta_2} is shape parameter one,
 \eqn{\beta_{0,1}} and \eqn{\beta_{0,2}} are the maximum response ,
 \eqn{\beta_3} is a scaling parameter two and  \eqn{\beta_4} is a shape parameter two.

 \item Quadratic model (\code{"qd"})
 \deqn{y=\beta_1 + \beta_2x + \beta_3x^2}
 where \eqn{\beta_1} is a constant, \eqn{\beta_2} is a linear coefficient
 and  \eqn{\beta_3} is the quadratic coefficient.

 \item Trapezium model (\code{"trapezium"})
 \deqn{y={\rm min}(\beta_1+\beta_2x, \beta_0, \beta_3 + \beta_4x)}
 where  \eqn{\beta_1} is the intercept one, \eqn{\beta_2} is the slope one,
 \eqn{\beta_0} is the maximum response, \eqn{\beta_3} is the intercept two
 and \eqn{\beta_3} is the slope two.

 \item Mitscherlich model (\code{"mit"})
 \deqn{y= \beta_0 - \beta_1*\beta_2^x}
 where \eqn{\beta_1} is the intercept, \eqn{\beta_2} is a shape parameter
 and \eqn{\beta_0} is the maximum response.

 \item Schmidt model (\code{"schmidt"})
 \deqn{y= \beta_0 + \beta_1(x-\beta_2)^2}
 where \eqn{\beta_1} is ascaling parameter, \eqn{\beta_2} is a
 shape parameter (x-value at maximum response ) and \eqn{\beta_0} is the
 maximum response .

 \item Custom model ("other")
 This option allows you to create your own model form using the function
 \code{function()}. The custom model should be assigned to the argument
 \code{equation}. Note that the parameters for the custom model should be
 \code{a} and \code{b} for a two parameter model; \code{a}, \code{b} and \code{c}
 for a three parameter model; \code{a}, \code{b}, \code{c} and \code{d} for a
 four parameter model and so on.
 }

The function \code{blbin()} utilities the optimization procedure of the
\code{optim()} function to determine the model parameters. There is a tendency
for optimization algorithms to settle at a local optimum. To remove the risk of
settling for local optimum parameters, it is advised that the function is run using
several starting values and the results with the smallest error (residue mean square)
can be taken as a representation of the global optimum.

The common errors encountered due to poor start values \enumerate{
\item function cannot be evaluated at initial parameters
\item initial value in 'vmmin' is not finite}
}
\examples{
x<-log(SoilP$P)
y<-SoilP$yield
start<-c(4,3,13.6, 35, -5)
bins<-c(1.6,4.74,0.314)

blbin(x,y, bins=bins, start=start,model = "trapezium", tau=0.99,
       xlab=expression("Phosphorus/ln(mg L"^-1*")"),
       ylab=expression("Yield/ t ha"^-1), pch=16,
       col="grey", bp_col="grey")

}
\references{
Casanova, D., Goudriaan, J., Bouma, J., & Epema, G. (1999). Yield gap analysis
in relation to soil properties in direct-seeded flooded rice.

Nelder, J.A. 1961. The fitting of a generalization of the logistic curve.
Biometrics 17: 89–110.

Phillips, B.F. & Campbell, N.A. 1968. A new method of fitting the von Bertelanffy
growth curve using data on the whelk. Dicathais, Growth 32: 317–329.

Schmidt, U., Thöni, H., & Kaupenjohann, M. (2000). Using a boundary line approach
to analyze N2O flux data from agricultural soils. Nutrient Cycling in Agro-ecosystems,
57, 119-129.
}
\author{
Chawezi Miti <chawezi.miti@nottingham.ac.uk>
}
