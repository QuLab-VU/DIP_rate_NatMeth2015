## **Code for Harris et al., Nature Methods 2016**
[doi:10.1038/nmeth.3852](http://dx.doi.org/10.1038/nmeth.3852)

All R code in this repository was written by _Darren Tyson_. All Python code was written by _Leonard Harris_.

Code for generating manuscript figures is provided in the `code_for_figs` directory. 
All experimental data needed to produce the figures is provided in the `data_for_figs` directory.

Standalone software for DIP rate estimation from experimental data is provided in the `dipDRC` directory.
An R script that runs an example application of the `dipDRC.r` code is provided in the `example_dipDRC` directory.

### **Instructions for running the R code**

The R statistical software package can be obtained free of charge from [www.R-project.org](http://www.R-project.org).

Required packages are:
```
deSolve
gplots
grid
drc
MASS
```
Assuming all required packages are installed, to generate all the graphs shown in the manuscript figures
and supplementary information, type the following in the R console:
```
source("[path_to_download_dir]/DIP_rate_NatMeth2016/code_for_figs/RcodeForFigs.r", chdir=TRUE)
```
where `[path_to_download_dir]` is the directory to which you downloaded this repository (e.g., `/Users/mycomp/git`).

An example application of the `dipDRC.r` code can be performed without downloading the entire Git repository
by copying all the code shown in:

([https://github.com/QuLab-VU/DIP_rate_NatMeth2016/blob/master/example_dipDRC/makeDRCexample.r](https://github.com/QuLab-VU/DIP_rate_NatMeth2016/blob/master/example_dipDRC/makeDRCexample.r))

#### *NOTE that running this code will download and install a number of R libraries automatically*

Alternatively, the code can be run from a local version using a similar format as for generting the manuscript figures:
```
source("[path_to_download_dir]/DIP_rate_NatMeth2016/example_dipDRC/makeDRCexample.r", chdir=TRUE)
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

### **Instructions for running the Python code**

The Python code for generating Supplementary Figure 7 ("Theoretical effects of variations around a mean cell seeding
density") is in the file `makeSuppFig7.py` within the `code_for_figs/Python` directory. Running the code requires installing `pysb`, 
a Python-based platform for biological modeling and simulation (see [www.pysb.org](http://pysb.org/)). PySB version 1.0.1
can be installed (Mac/Linux) by opening a command prompt and typing
```
sudo pip install pysb
``` 
PySB also requires installing BioNetGen ([www.bionetgen.org](http://www.bionetgen.org)). Please see 
[docs.pysb.org/en/latest/installation.html](http://docs.pysb.org/en/latest/installation.html) for detailed installation 
instructions.

Once PySB and all its dependencies are installed, `makeSuppFig7.py` can be run by simply typing
```
python makeSuppFig7.py
```
The output is three figures: `timecourses.pdf` (Supplementary Fig. 7a), `boxplots.pdf` (Supplementary Fig. 7b), and `distributions.pdf`
(Supplementary Fig. 7c,d). Functionally, the code runs 10^5 control simulations and 10^4 drug-treated simulations for each
of two scenarios: (i) variations in cell seeding density alone, and (ii) variations in both sampling time and cell seeding
density. This amounts to a total of 2\*10^5+2\*10^4=220,000 simulations (i.e., it may take a while to finish).
