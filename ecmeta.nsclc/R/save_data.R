#' @export
save_data <- function() {
  psmodel_specification <- read.csv("data-raw/psmodel-specification.csv")
  save(psmodel_specification, file = "data/psmodel_specification.rda", compress = "bzip2")
}