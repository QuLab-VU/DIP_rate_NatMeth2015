##**Code for generating figures for Harris et al., Nature Methods 2016**

All code in this repository was written by _Darren Tyson_ (R) and _Leonard Harris_ (Python).

The R statistical software package is required to produce figures and can be obtained for free from [www.R-project.org](http://www.R-project.org)

Also included in this repository is all necessary data in comma-separated value (.csv) format.

To generate all figures, R must be in the directory in which this README file is found.

Assuming you have all the required libraries installed, at the R command line prompt simply type:
```
source("./code_for_figs/R/codeForFigs.R",chdir=TRUE)
```

To see an example of how the `dipDRC` function takes raw cell count data, 
estimates DIP rates and produces dose-response curves, at the R command line prompt
type:
```
source("./example_dipDRC/makeDRCexample.R",chdir=TRUE)
```

The data used in this example is in the file named: `dipDRC_example_data.csv`
It is the raw breast cancer cell line (cell counts) data used to make Figure 4: MDA-MB-231 cells treated with 
rotenone or phenformin. The data are from a single experiment with two technical replicates and are 
structured as a seven-column matrix with the following column headers:
1) `time`
2) `cell.count`
3) `cell.line`
4) `drug`
5) `conc`
6) `well`
7) `expt.date`

The output of the `dipDRC` function is a graphics window with two plots, one for each cell line and drug condition
and a list of the two `drm` (dose-response model) objects, described in detail in the package information of the
`drc` library: [https://cran.r-project.org/web/packages/drc/drc.pdf](https://cran.r-project.org/web/packages/drc/drc.pdf)
