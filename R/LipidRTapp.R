# LipidRTapp()
#' LipidRT shiny app
#'
#' Interactive UI for filtering and RT prediction LipidRT
#'
#' @details To provide a more interactive and user-friendly experience, a Shiny
#' application was developed to implement the workflow through a graphical interface
#' in R. This application enables users to execute all steps of the analysis in a
#' visual and intuitive manner.
#'
#' @examples
#' \dontrun{
#' # example data files can be download from github.com/maialba3/LipidMSv2.0_exampleFiles
#'
#' library(LipidRT)
#' LipidRTapp()
#' }
#'
#' @author Celia Regal Marques <celia_regal@iislafe.es>
LipidRTapp <- function() {
  for (pkg in c("shiny", "shinythemes", "shinyjs", "plotly", "zip")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' must be installed to use this function."),
           call. = FALSE)
    }
  }
  if (interactive()){
    shiny::runApp(system.file('LipidRTapp', package='LipidRT'))
  }else {
    stop("Only available for interactive sessions")
  }
}
