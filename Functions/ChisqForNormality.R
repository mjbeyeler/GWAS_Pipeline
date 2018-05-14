ChisqForNormality <- function(x, nbins=0) {
  # This function calculates a chi-squared goodness-of-fit
  # * Plots actual and theoretical densities.
  # * Returns a p-value, estimating the probability that the source data are normally distributed
  
  # By default, it uses the Freedman-Diaconis method to determine the 'optimal' number of bins. but you can also enter a custom number you want.
  
  # Usually, a goodness of fit evaluation loses 1 degree of freedom because the number of total observations must be the same.
  # Here, we lose 2 additional df for the 2 parameters we have to estimate from the data:
  # * the mean
  # * the standard deviation
  # which we use to generate the ideal normal distribution.
  
  x.length <- length(x)
  x.min <- min(x)
  x.max <- max(x)
  x.mean <- mean(x)
  x.sd <- sd(x)
  
  if (nbins!=0)
    my.breaks <- seq(x.min, x.max, length=nbins+1)
  else
    my.breaks <- 'FD'
  
  h=hist(x, breaks=my.breaks, plot=F)
  break.width <- h$mids[2]-h$mids[1]
  
  theoretical.density <- dnorm(h$mids, x.mean, x.sd)
  theoretical.pdf <- theoretical.density * break.width
  pdf.missing <- 1 - sum(theoretical.pdf)
  theoretical.pdf <- c(pdf.missing/2, theoretical.pdf, pdf.missing/2)
  theoretical.counts <- theoretical.pdf * x.length
  counts <- c(0, h$counts, 0)
  
  # print(counts)
  # print(theoretical.counts)
  
  
  # Before counting the X², make sure that each bin contains at least 5 samples.
  # The following loop combines bins until all of them contain at least 
  while(min(counts)<5 | min(theoretical.counts)<5) {
    for (i in 1:length(counts)) {
      if (i < length(counts)) {
        if(counts[i] < 5 | theoretical.counts[i] < 5) {
          counts[i+1] <- counts[i+1] + counts[i]
          theoretical.counts[i+1] <- theoretical.counts[i+1] + theoretical.counts[i]
          counts <- counts[-i]
          theoretical.counts <- theoretical.counts[-i]
          break
        }
      }
      # This part is needed in case the last bin of the histocram contains less than 5 counts.
      else {
        if(counts[i] < 5 | theoretical.counts[i] < 5) {
          counts[i-1] <- counts[i-1] + counts[i]
          theoretical.counts[i-1] <- theoretical.counts[i-1] + theoretical.counts[i]
          counts <- counts[-i]
          theoretical.counts <- theoretical.counts[-i]
          break
        }
      }
    }
  }
  plot(counts, type='l')
  lines(theoretical.counts, col=2)
  # Doing the actuall chi-squared goodness-of-fit test:
  X2 <- sum((counts-theoretical.counts)^2 / theoretical.counts)
  # For normal distribution:
  bins.final <- length(counts)
  df <- bins.final - 3
  print(paste('X2 =', X2, 'bins =', bins.final, 'df =', df))
  return(pchisq(X2, df, lower.tail=F))
}