if(!require("knitr", character.only=T))
  install.packages("knitr", repos="http://cran.us.r-project.org")

purl(commandArgs(trailingOnly=TRUE))
