# Running this file replicates and retrieves all results reported in the paper
# Since the analysis takes a long time to run, two files are created to monitor progress:
# (1) ecmeta.nsclc/vignettes/psweight_out reports progress when generating propensity score weights 
# for each study (the genetic matching algorithm is slow)
# (2) simulation/output/sim-out tracks progress while simulating the data and evaluating
# the model(s) with Monte Carlo replications

run_simulation <- function() {
  rmarkdown::render(input = "simulation/simulations.Rmd", output_dir = "simulation")
}

run_nsclc <- function() {
  pkgdown::build_site("ecmeta.nsclc")
}

write_txtstats <- function(txt) {
  # Convert statistics to data frame
  txtstats <- data.frame(do.call(rbind, txt))
  
  # Output to text file to input into latex
  txtstats$def <-  "\\def"
  names(txtstats)[1] <- "value"
  txtstats$value <- as.character(txtstats$value)
  txtstats <- data.frame(def = txtstats$def, name = rownames(txtstats), value =  txtstats$value)
  txtstats$output <- paste(txtstats[, 1], " ", "\\", txtstats[, 2],
                           "{", txtstats[, 3], "}", sep = "")
  fileConn <- file("txtstats.txt")
  writeLines(txtstats$output, fileConn)
  close(fileConn)
}

get_results <- function() {
  # aNSCLC example
  ## Figures
  fig_names <- file.path(
    "ecmeta.nsclc/vignettes/figs",
    list.files("ecmeta.nsclc/vignettes/figs")
  )
  file.copy(from = fig_names, to = "figs/", overwrite = TRUE)
  
  ## Tables
  tbl_names <- file.path(
    "ecmeta.nsclc/vignettes/tables",
    list.files("ecmeta.nsclc/vignettes/tables")
  )
  file.copy(from = tbl_names, to = "tables/", overwrite = TRUE)
  
  ## Inline text
  write_txtstats(readRDS("ecmeta.nsclc/vignettes/txtstats.rds"))
  
  # Simulation
  ## Figures
  fig_names <- file.path(
    "simulation/figs",
    list.files("simulation/figs")
  )
  file.copy(from = fig_names, to = "figs/", overwrite = TRUE)
  
  
  invisible()
}

replicate <- function(what = c("example", "simulation"), get = TRUE) {
  
  what <- match.arg(what, several.ok = TRUE)
  cat("Running replication script: \n", file = "replication-out.txt")
  
  # Simulation
  if ("simulation" %in% what) {
    cat("Starting simulation: \n", file = "replication-out.txt", append = TRUE)
    rmarkdown::render(input = "simulation/simulations.Rmd", output_dir = "simulation")
    cat("Completed simulation: \n", file = "replication-out.txt", append = TRUE)
  }

  # aNSCLC example
  if("ecdata" %in% rownames(installed.packages())) {
    if ("example" %in% what) {
      cat("Starting aNSCLC analysis: \n", file = "replication-out.txt", append = TRUE)
      setwd("ecmeta.nsclc")
      pkgdown::build_site()
      setwd("..")
      cat("Completed aNSCLC analysis: \n", file = "replication-out.txt", append = TRUE)
    } 
  } else{
    message(paste0(
      "The aNSCLC example was not replicated because the Flatiron Health ",
      "database required to build the external control cohorts is not ",  
      "available. If you would like to build the external control cohorts, ",
      "you will need to install the Roche package 'ecdata'."
    ))
  }
  
  # Get results
  if (get) get_results()
}

replicate(what = c("simulation", "example"), get = TRUE)