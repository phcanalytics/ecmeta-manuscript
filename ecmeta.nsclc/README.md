# Decision Making and Error Control for External Control Analyses in NSCLC
This `R` package contains code for a proof of concept study in advanced non-small cell lung cancer (NSCLC) of a meta-analytic framework to adjust external control studies for additional bias and variability due to their non-randomized design. It is not meant to be a general purpose library; instead, it simply aims to facilitate transparent, reproducible, and modular analysis. The analyses estimate the parameters of the meta-analytic model and illustrate use for a new hypothetical single-arm trial. 

For a general purpose package, see [`ecmeta`](https://pages.github.roche.com/RWDScodeshare/ecmeta/), which describes, assesses, and implements the framework for general use cases.

## Methodology
Details of the meta-analytic framework can be found [here](https://pages.github.roche.com/RWDScodeshare/ecmeta/articles/methodology.html).

## Illustration
The framework can only be applied to a new trial after estimation of parameters using relevant historical data from RCTs and external controls. This project estimates parameters relevant for a new trial evaluating treatment of patients with advanced NSCLC and illustrates application of the framework to this new trial. The prespecified statistical analysis plan for parameter estimation can be found [here](https://pages.github.roche.com/RWDScodeshare/ecmeta.nsclc/sap/sap.pdf). 14 pairwise comparisons from 11 separate clinical trials are used in the analysis.

## Intallation
The `ecmeta.nsclc` package can be installed from the `R` console with:

```{r}
system("R CMD INSTALL ecmeta.nsclc")
```
