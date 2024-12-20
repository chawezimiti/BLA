% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Perimeter_Area_calc.R
\name{AP}
\alias{AP}
\title{Calculate perimeter and area of data scatter}
\usage{
AP(points)
}
\arguments{
\item{points}{A dataframe with two columns.}
}
\value{
The area and perimeter covered bivariate data points, and the selected
  boundary points
}
\description{
This function determines the perimeter and area around the boundary points
}
\examples{
x<-data.frame(x=evapotranspiration$`ET(mm)`,y=evapotranspiration$`yield(t/ha)`)
AP(x)

}
\keyword{internal}
