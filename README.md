##**Code for generating figures for Harris et al., Nature Methods 2016**

All R code in this repository was written by _Darren Tyson_. All Python code was written by _Leonard Harris_.

Code for generating manuscript figures is provided in the `code_for_figs` directory. 
All experimental data needed to produce the figures is provided in the `data_for_figs` directory.

Standalone software for DIP rate estimation from experimental data is provided in the `dipDRC` directory.
An R script that runs an example application of the `dipDRC.r` code is provided in the `example_dipDRC` directory.

=== Instructions for running the R code ===

The R statistical software package can be obtained free of charge from [www.R-project.org](http://www.R-project.org)

Required packages are
```
deSolve
gplots
grid
drc
MASS
```
Assuming all required packages are installed, in the R console type:
```
setwd("[path_to_download_dir]/DIP_rate_NatMeth2016/code_for_figs/R")
source("codeForFigs.R")
```
where `[path_to_download_dir]` is the directory to which you downloaded this repository.

Similarly, to run the example using the `dipDRC.r` code that takes raw cell count data,  
estimates DIP rates, and produces dose-response curves, type:
```
setwd("[path_to_download_dir]/DIP_rate_NatMeth2016/example_dipDRC")
source("makeDRCexample.R")
```
The data used in this example are the raw breast cancer (MDA-MB-231) cell counts for rotenone
and phenformin treatments. They can be found in the file `example_dipDRC/dipDRC_example_data.csv`. The data 
are from a single experiment with two technical replicates. The data file is structured as a 
seven-column matrix with the following headers:
1) `time`
2) `cell.count`
3) `cell.line`
4) `drug`
5) `conc`
6) `well`
7) `expt.date`

The output is a graphics window with two plots, one for each cell line and drug condition,
and a list of the two `drm` (dose-response model) objects, described in detail in the package information of the
`drc` library [https://cran.r-project.org/web/packages/drc/drc.pdf](https://cran.r-project.org/web/packages/drc/drc.pdf)

=== Instructions for running the Python code ===

_Coming soon_
