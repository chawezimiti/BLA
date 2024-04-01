#' Print class wm
#'
#' Prints an output with 5 elements.
#'
#' @param x Print object.
#' @param ... Other parameters associated with the \code{print()} function.
#'
#' @keywords internal
#' @export
#'
print.wm <- function(x,...) {
  print(x[c(1:4)])}  # define a new print function for the cm class
