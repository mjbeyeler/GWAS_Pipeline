s <- 43

source('Helper_Scripts/Reprobucibility_Functions.R')
.BuildReproducibleEnvironment('2018-01-01')

.LIST.OF.PACKAGES <- c(
  'data.table',           #
  'tictoc',               # 
  'lme4',                 # 
  'reticulate',           # used to properly being able to switch between python and r in R Markdown
  'icesTAF'               # dos2unix function
)
.LoadPackages(.LIST.OF.PACKAGES)
s2 <- 44
sessionInfo()