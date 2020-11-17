plot.snp.by.function <- function(all.gene.info, pareto.col = "#9ecae1", max.fp = 0.2){
	
	emma <- as.numeric(all.gene.info[,"EMMA.p"])
	fp <- as.numeric(all.gene.info[,"FP"])

	pref <- high(emma) * high(fp)

	plot(emma, fp)
	plot_front(data.frame(cbind(emma, fp)), pref)

	
	plot.new()
	plot.window(xlim = c(min(emma), max(emma)), ylim = c(min(fp), max(fp)))
	plot_front(data.frame(cbind(emma, fp)), pref, col = pareto.col, lwd = 3)

	if(!is.null(max.fp)){
		abline(h = -log10(max.fp), lty = 2, col = "darkgray")	
	}
	par(xpd = TRUE)
	text(emma, fp, labels = rownames(all.gene.info), cex = 0.7)
	par(xpd = FALSE)
	axis(1);axis(2)
	mtext("-log EMMA p value", side = 1, line = 2)	
	mtext("-log SVM FP rate", side = 2, line = 2)
				
}