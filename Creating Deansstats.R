
# First make sure you have a clean working environment

rm(list = ls())



# Install and load devtools

install.packages('devtools')
library(devtools)

# Set working directory inside the folder 'Deansstats'. 
#Then run this function to create the R-package 'Deansstats'

devtools::build()


# Set working directory outside the folder 'Deansstats"  to the folder
# where 'Deansstats' is located. Then install the package and other package
# dependencies
install.packages("MASS")
install.packages("Matrix")
install.packages("ggplot2")
devtools::install("Deansstats")

# Now you can load all packages and use all functions from Deansstats. 

library(MASS)
library(Matrix)
library(ggplot2)
library(Deansstats)

# Try to run the help function for one of the functions to ensure the package
# is properly loaded

?RSA

# If it works, great! You are all set to work with the package. 
# If it does not work, close Rstudio. Then, reload the package and try again.

library(Deansstats)
?RSA

#Now, it should work.



