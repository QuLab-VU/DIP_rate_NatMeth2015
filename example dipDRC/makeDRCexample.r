# if devtools library not insalled, uncomment next line
# install.packages('devtools')

# libraries necessary to load source files directly from GitHub
library(devtools)
library(RCurl)

# pull source code directly from GitHub
SourceURL <- "https://raw.github.com/QuLab-VU/DIP_rate_NatMeth2016/master/dipDRC/dipDRC.r"
source_url(SourceURL)

# download example data to local directory 
# comment out 2 lines below to load your own data
download.file("https://raw.githubusercontent.com/QuLab-VU/DIP_rate_NatMeth2016/master/dipDRC/dipDRC_example_data.csv", 
    destfile = "dipDRC_example_data.csv", method = "curl")
   

# load example data into object 'd'
# if you want to load your own data, 
# replace 'dipDRC_example_data.csv' with file name of your data
d <- read.csv('dipDRC_example_data.csv')


# identify DIP rates for each cell line and condition and produce a dose-response
# curve for each. Model fits to data are stored into the variable 'dose.resp'
dose.resp	<-	dipDRC(d, type='confidence')





