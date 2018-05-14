NormalityHistogram <- function(x) {
  x.length <- length(x)
  x.min <- min(x)
  x.max <- max(x)
  x.mean <- mean(x)
  x.sd <- sd(x)
  
  nbins <- 25
  my.breaks <- seq(x.min, x.max, length=nbins+1)
  break.width <- my.breaks[2]-my.breaks[1]
  

  # First plot without density:
  # The prob=T part is necessary for the densiity curve to fit the histogram.
  h=hist(x, breaks=my.breaks, prob=T)
  
  # For testing
  #print(h)

  # Second plot, with density and normal approximation
  hist(x, breaks=my.breaks, prob=T)
  lines(density(x), col=2)
  
  # Approximating data as normal distribution from the midpoints of the histogram
  x.as.normal <- dnorm(h$mids, x.mean, x.sd)
  lines(h$mids, x.as.normal, col=3)
  # print('wha')
}