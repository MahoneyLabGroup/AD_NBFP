smooth.inv.cdf <- function(data,samp.freq=1000,min.val=NULL,max.val=NULL){

# Load libraries
library(splines)

# Set boundaries for the smoothing
if(is.null(min.val)){
	min.val = min(data) - 1/samp.freq
	}

if(is.null(max.val)){
	max.val = max(data) + 1/samp.freq
	}

# Smooth the data
den = density(data,n=samp.freq,from=min.val,to=max.val)

# Get cumulative density function (CDF)
den.cdf = cumsum(den$y)/sum(den$y)

# Get interpolant for inverse CDF
inv.cdf.fit = splinefun(den.cdf,den$x)

return(inv.cdf.fit)

}