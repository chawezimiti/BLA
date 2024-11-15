#' Calculate perimeter and area of data scatter
#'
#' This function determines the perimeter and area around the boundary points
#'
#' @param points A dataframe with two columns.
#'
#' @returns The area and perimeter covered bivariate data points, and the selected
#'   boundary points
#' @keywords internal
#' @export
#'
#' @examples
#' x<-data.frame(x=evapotranspiration$`ET(mm)`,y=evapotranspiration$`yield(t/ha)`)
#' AP(x)
#'
AP<- function(points) {

  perimeter <- function(points) {
    # Initialize perimeter
    perimeter <- 0

    # Loop over pairs of consecutive points and calculate Euclidean distance
    for (i in 1:(nrow(points) - 1)) {
      dx <- points[,1][i+1] - points[,1][i]
      dy <- points[,2][i+1] - points[,2][i]
      perimeter <- perimeter + sqrt(dx^2 + dy^2)
    }

    # Add distance between the last and first point to close the polygon
    dx <- points[,1][1] - points[,1][nrow(points)]
    dy <- points[,2][1] - points[,2][nrow(points)]
    perimeter <- perimeter + sqrt(dx^2 + dy^2)

    return(perimeter)
  }

  # Calculate the perimeter
  perimeter <- perimeter(points) #chage back to point

  # Calculate the area of the polygon using the Shoelace formula

  area <- function(points) {
    # hull_points: A matrix or data frame with two columns (x and y coordinates)

    n <- nrow(points)  # Number of points in the hull
    area <- 0

    for (i in 1:(n - 1)) {
      area <- area + (points[i, 1] * points[i + 1, 2] - points[i + 1, 1] * points[i, 2])
    }

    # Add the last term (closing the polygon)
    area <- area + (points[n, 1] * points[1, 2] - points[1, 1] * points[n, 2])

    # Take the absolute value and divide by 2
    area <- abs(area) / 2
    return(area)
  }

  # Calculate the area

  area <- area(points)


  return(list(Perimeter=perimeter,Area=area))
}

