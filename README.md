##**Code for Harris et al., Nature Methods 2016**

All R code in this repository was written by _Darren Tyson_. All Python code was written by _Leonard Harris_.

Code for generating manuscript figures is provided in the `code_for_figs` directory. 
All experimental data needed to produce the figures is provided in the `data_for_figs` directory.

Standalone software for DIP rate estimation from experimental data is provided in the `dipDRC` directory.
An R script that runs an example application of the `dipDRC.r` code is provided in the `example_dipDRC` directory.

###**Instructions for running the R code**

The R statistical software package can be obtained free of charge from [www.R-project.org](http://www.R-project.org).

Required packages are:
```
deSolve
gplots
grid
drc
MASS
```
Assuming all required packages are installed, to run the code type the following in the R console:
```
setwd("[path_to_download_dir]/DIP_rate_NatMeth2016/code_for_figs/R")
source("codeForFigs.r")
```
where `[path_to_download_dir]` is the directory to which you downloaded this repository (e.g., `/Users/mycomp/git`).

Similarly, to run the example application of the `dipDRC.r` code type:
```
setwd("[path_to_download_dir]/DIP_rate_NatMeth2016/example_dipDRC")
source("makeDRCexample.r")
```
The output is a graphics window with two plots, one for each cell line and drug condition, and a list of the two 
`drm` (dose-response model) objects, described in detail in the package information for the `drc` library 
([https://cran.r-project.org/web/packages/drc/drc.pdf](https://cran.r-project.org/web/packages/drc/drc.pdf)).

The data used in this example are the raw breast cancer (MDA-MB-231) cell counts for rotenone
and phenformin treatments. They can be found in the file `dipDRC_example_data.csv` within the `example_dipDRC` directory. 
The data are from a single experiment with two technical replicates. The data file is structured as a 
seven-column matrix with the following headers:
1) `time`
2) `cell.count`
3) `cell.line`
4) `drug`
5) `conc`
6) `well`
7) `expt.date`

###**Instructions for running the Python code**

The Python code for generating Supplementary Figure 7 ("Theoretical effects of variations around a mean cell seeding
density") is in the file `makeSuppFig7.py` within the `code_for_figs/Python` directory. Running the code requires installing `pysb`, 
a Python-based platform for biological modeling and simulation (see [www.pysb.org](http://pysb.org/)). This can be done
by opening a command prompt and typing
```
sudo pip install pysb
``` 
(this will install pysb version 1.0.1). `makeSuppFig7.py` can then be run by simply typing
```
python makeSuppFig7.py
```
The output is three figures: `timecourses.pdf` (Supplementary Fig. 7a), `boxplots.pdf` (Supplementary Fig. 7b), and `distributions.pdf`
(Supplementary Fig. 7c,d). Functionally, the code runs 10^5 control simulations and 10^4 drug-treated simulations for each
of two scenarios: (i) variations in cell seeding density alone, and (ii) variations in both sampling time and cell seeding
density, for a total of 2*10^5+2*10^4=220,000 simulations (i.e., it may take a while to finish).