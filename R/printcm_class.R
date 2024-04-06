#' Print class cm
#'
#' This is an S3 print function that Prints only the first 4 elements of
#' cm class objects.
#'
#' @param x Print object.
#' @param ... Other parameters associated with the \code{print()} function.
#' @returns A object containing only the first four items.
#' @keywords internal
#' @export
#' @examples
#' numbers<- 1:10
#' class(numbers)<-"cm"
#' numbers
#'
print.cm <- function(x,...) {
  print(x[c(1:4)]) } # define a new print function for the cm class
