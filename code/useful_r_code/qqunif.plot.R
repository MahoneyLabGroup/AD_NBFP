
qqunif.plot<-function(pvalues, plot.label = "P Value qq plot"){

	xlab=expression(paste("Expected (",-log[10], " p-value)"))
	ylab=expression(paste("Observed (",-log[10], " p-value)")) 
	
			
	exp.dist <- runif(length(pvalues))
	log.exp <- -log10(exp.dist)
	log.obs <- -log10(pvalues)

	qqplot(log.exp, log.obs, xlab = xlab, ylab = ylab, main = plot.label)
	abline(0,1)

}