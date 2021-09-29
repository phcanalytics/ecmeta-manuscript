#' @export
run_app <- function() {
  f <- system.file("app", "ps", package = "ecmeta.nsclc")
  if (f == "") {
    stop("Could not find Shiny app directory. Try re-installing `ecmeta.nsclc`.", 
         call. = FALSE)
  }
  shiny::runApp(f, display.mode = "normal")
}