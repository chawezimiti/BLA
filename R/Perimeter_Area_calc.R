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

  #Combine x and y into a data frame

  data <- data.frame(x = points[,1], y = points[,2])

  #Separate data frames for max y and min y for each x
  max_data <- data[order(data$x, -data$y), ]  # Sort by x ascending, y descending
  max_data <- max_data[!duplicated(max_data$x), ]  # Keep only max y for each x
  max_data$selection <- "max"  # Add a column to indicate max selection

  min_data <- data[order(data$x, data$y), ]  # Sort by x ascending, y ascending
  min_data <- min_data[!duplicated(min_data$x), ]  # Keep only min y for each x
  min_data$selection <- "min"  # Add a column to indicate min selection

  #Combine max and min data into a single data frame

  result <- rbind(max_data, min_data)
  point<- result

  #Create the polygon by ordering max points (in ascending x) and min points (in descending x)
  polygon_points <- rbind(
    result[result$selection == "max", ],  # Max points in ascending x
    result[result$selection == "min", ][nrow(result[result$selection == "min", ]):1, ]  # Min points in descending x
    )

  # Calculate the perimeter of the polygon

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
  perimeter <- perimeter(point)

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

