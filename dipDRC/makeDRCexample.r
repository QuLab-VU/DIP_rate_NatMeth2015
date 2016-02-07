#install.packages('devtools')
# libraries necessary to load source files directly from GitHub
library(devtools)
library(RCurl)

# pull source code directly from GitHub
SourceURL <- "https://raw.github.com/QuLab-VU/DIP_rate_NatMeth2016/master/dipDRC/dipDRC.r"
source_url(SourceURL)

# Download example data (comment out if loading your own data)
# make temporary file for data to be downloaded
temp <- tempfile()
# download example data to temp file 
download.file("https://raw.githubusercontent.com/QuLab-VU/DIP_rate_NatMeth2016/master/dipDRC/dipDRC_example_data.csv", 
    destfile = temp, method = "curl")
   
# load example data into object 'd'
# if you want to load your own data, 
# replace temp with file name of your data
d <- read.csv(temp)

# if using example data, remove temporary file
unlink(temp)
rm(temp)

# identify DIP rates for each cell line and condition and produce a dose-response
# curve for each. Model fits to data are stored into the variable 'dose.resp'
dose.resp	<-	dipDRC(d, type='confidence')





