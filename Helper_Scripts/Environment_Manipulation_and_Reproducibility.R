if(!require('checkpoint', character.only=T))
  install.packages('checkpoint', repos="http://cran.us.r-project.org")


# LOADING PACKAGES FUNCTION

# Function that loads list of packages, and installs them if necessary.

.load_packages <- function(package.list) {
  invisible(lapply(package.list, function(not.a.package) {if (!require(not.a.package, character.only=T)) {
    # if(not.a.package %in% loadedNamespaces()) {
    #   unloadNamespace(not.a.package)
    #   print(not.a.package)
    # }
    install.packages(not.a.package, repos="http://cran.us.r-project.org");
    if(!require(not.a.package, character.only=T))
      stop(paste(not.a.package, "Package not found"))}}
  ))
}



# REMOVING ALL EXISTING CHECKPOINTS

.remove_all_checkpoints <- function() {
  sapply(checkpoint::checkpointArchives(), checkpoint::checkpointRemove)
}



# BUILD A REPRODUCIBLE ENVIRONMENT

.build_reproducible_environment <- function(PROJECT.SNAPSHOT.DATE = substr(R.version.string, 18, 27),
                                          PROJECT.VERSION       = paste(R.version$major, R.version$minor, sep="."),
                                          SCAN.FOR.PACKAGES     = TRUE) {
  
  current.version <- paste(R.version$major, R.version$minor, sep=".")
  if(current.version != PROJECT.VERSION)
    cat(paste("Warning: Your R version (", current.version,
              ") differs from the version this project was programmed in (", PROJECT.VERSION, "),\n",
              "which might be causing inconveniences.\n\n",
              sep=""))
  
  
  # SETTING UP A CHECKPOINT IF NECESSARY
  # A new checkpoint is set if the user
  # (a) is using CRAN, or
  # (b) is using a idfferent version of MRAN (and thus a different snapshot date).
  
  tmp <- checkpoint::setSnapshot()[1]
  current.project.date <- substr(tmp, nchar(tmp) - 9, nchar(tmp))
  # Note: If current.project.date will be a date. If using CRAN, and having no chekckpoint loaded, it will be something like
  # "tudio.com/" or something similar.
  
  # If the current checkpoint differs from the one to be used for a project - or if no checkpoint is
  # loaded at all - it will be changed to the desired date, and installing all necessary packages for that snapshot date.
  if(current.project.date != PROJECT.SNAPSHOT.DATE) {
    if (!dir.exists(paste("~/.checkpoint/", PROJECT.SNAPSHOT.DATE, sep=""))) {
      if (!dir.exists("~/.checkpoint/"))
        dir.create("~/.checkpoint/")
      dir.create(paste("~/.checkpoint/", PROJECT.SNAPSHOT.DATE, sep=""))
    }
    setwd("Helper_Scripts/Empty_Folder")
    checkpoint::checkpoint(PROJECT.SNAPSHOT.DATE, scanForPackages= SCAN.FOR.PACKAGES)
    checkpoint::setSnapshot(PROJECT.SNAPSHOT.DATE)
    setwd("../..")
    }
  
  # To speed up matrix calculations through multithreading for MRAN users
  if(exists("Revo.version")) {
    library("RevoUtilsMath")
    library("RevoUtils")
  }
}



# CLEARING ENVIRONMENT FUNCTION

# This function is as close as I got to reproduce what the Ctrl+Shift+F10 command does (which starts a fresh RStudio session).
# I suspect some of this code might do some unwanted things, thus I abandon this function for now.
# -> Use Ctrl+Shift+F10 (or Session -> Restart R) instead, when wanting to start from scratch
#    This is also what hadley recommends: https://community.rstudio.com/t/first-line-of-every-r-script/799/36?u=mbey

# Even more so: Knitting apparently can't be disturbed by the environment, as it runs in a sandboxed area. So, cleaning the environment isn't important. Only clean in cases of emergency.

# .clear_environment <- function() {
#   
#   # Detatches and unloads all the attached packages except the base ones,
#   # which should not be detached anyway
#   basic.packages <- c("package:stats", "package:graphics", "package:grDevices",
#                       "package:utils", "package:datasets", "package:methods",
#                       "package:base")
#   package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#   package.list <- setdiff(package.list,basic.packages)
#   if (length(package.list) > 0) 
#     for (package in package.list) {
#       detach(package, character.only=TRUE, unload=TRUE)
#       print(paste("Package ", package, " detached.", sep = ""))
#     }
#   
#   # Cleans the global environment,
#   rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
#   
#   # Clears all displayed plots on the standard plot device, if there are any
#   while(!is.null(dev.list())) {
#     dev.off(dev.list()[1])}
# }