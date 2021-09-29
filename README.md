# A meta-analytic framework to adjust for bias in external control studies
This repository contains code for a paper that introduces a meta-analytic approach for adjusting external control studies for additional bias and variability due to their non-randomized design. Code for the advanced non-small cell lung cancer (aNSCLC) example is available in the `ecmeta.nsclc` directory and a corresponding website with results is available [here](https://pages.github.roche.com/incertid/ecmeta-manuscript/ecmeta.nsclc/docs/). Code for the simulation study is available in the `simulation` directory and output from the simulations is available [here](https://pages.github.roche.com/incertid/ecmeta-manuscript/simulation/simulations.html).

## Replication
The entire analysis can be replicated by running the file `replicate.R`, which runs code to replicate the simulation study and the aNSCLC example, and retrieves all results (figures, tables, and inline text) that are reported in the manuscript.

## Data availability 
The aNSCLC analyses can only be reproduced with access to the Flatiron Health database, since it was used to construct the external control cohorts. Unfortunately, we are not permitted to share this data.

## Dependencies
`R` package dependencies are managed through the [`renv`](https://rstudio.github.io/renv/articles/renv.html) package. You can view all packages and their versions in the lockfile [`renv.lock`](renv.lock). All required packages and the appropriate versions can be installed with `renv::restore()`. 

Note that an `R` package was created for the aNSCLC example, which can be installed with:

```{r}
system("R CMD INSTALL --build  --no-multiarch --with-keep.source ecmeta.nsclc")
```



