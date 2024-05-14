#' Check and Load Tcl/Tk
#'
#' This function checks if the Tcl/Tk environment is available and loads the `tcltk`
#' package if possible.
#' It returns TRUE if `tcltk` is successfully loaded, otherwise FALSE.
#' @return Logical value indicating whether `tcltk` was successfully loaded.
#' @keywords internal
#' @import tcltk
#' @export
#'
check_and_load_tcltk <- function() {
  if (interactive() || Sys.getenv("DISPLAY") != "") {
    if (requireNamespace("tcltk", quietly = TRUE)) {
      library(tcltk)
      return(TRUE)
    } else {
      message("The 'tcltk' package is not available.")
      return(FALSE)
    }
  } else {
    message("DISPLAY variable is not set, 'tcltk' functionality is unavailable.")
    return(FALSE)
  }
}
