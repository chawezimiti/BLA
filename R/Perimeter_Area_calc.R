#' Calculate perimeter and area of data scatter
#'
#' This function determines the perimeter and area around the boundary points
#'
#' @param dat A dataframe with two columns.
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
AP<- function(dat) {

  # Combine x and y into a data frame
  data <- data.frame(x = dat[,1], y = dat[,2])

  # Separate data frames for max y and min y for each x
  max_data <- data[order(data$x, -data$y), ]  # Sort by x ascending, y descending
  max_data <- max_data[!duplicated(max_data$x), ]  # Keep only max y for each x
  max_data$selection <- "max"  # Add a column to indicate max selection

  min_data <- data[order(data$x, data$y), ]  # Sort by x ascending, y ascending
  min_data <- min_data[!duplicated(min_data$x), ]  # Keep only min y for each x
  min_data$selection <- "min"  # Add a column to indicate min selection

  # Combine max and min data into a single data frame
  result <- rbind(max_data, min_data)

  # Create the polygon by ordering max points (in ascending x) and min points (in descending x)
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
      dx <- points$x[i+1] - points$x[i]
      dy <- points$y[i+1] - points$y[i]
      perimeter <- perimeter + sqrt(dx^2 + dy^2)
    }

    # Add distance between the last and first point to close the polygon
    dx <- points$x[1] - points$x[nrow(points)]
    dy <- points$y[1] - points$y[nrow(points)]
    perimeter <- perimeter + sqrt(dx^2 + dy^2)

    return(perimeter)
  }

  # Calculate the area of the polygon using the Shoelace formula
  area <- function(points) {
    n <- nrow(points)
    area <- 0
    for (i in 1:(n - 1)) {
      area <- area + (points$x[i] * points$y[i + 1] - points$x[i + 1] * points$y[i])
    }
    # Add the last term to close the polygon
    area <- area + (points$x[n] * points$y[1] - points$x[1] * points$y[n])

    return(abs(area) / 2)
  }

  # Calculate the perimeter and area
  perimeter <- perimeter(result)
  area <- area(result)


  return(list(Perimeter=perimeter,Area=area,Polygon=polygon_points))
}

