plot.final.results <- function(all.gene.info, snps, genetic.regions, sig.p = 0.05, max.fp = 0.2, value = "gene.final.score", label = "Final Gene Score"){
	
	library(RColorBrewer)
	region.cols <- brewer.pal(4, "Accent")
	
	layout(matrix(c(1,2), ncol = 1))
	par(mar = c(0, 4,4,2))

	snp.pos <- c(snps[,"Position"])
	all.nlp <- -log10(snps[,"p.value"])
	x.lim <- c(min(snp.pos, na.rm = TRUE), max(snp.pos, na.rm = TRUE))
	y.lim <- c(min(all.nlp), max(all.nlp))

	plot.new()
	plot.window(xlim = x.lim, ylim = y.lim)
	for(i in 1:nrow(genetic.regions)){
		draw.rectangle(genetic.regions[i,1], genetic.regions[i,2], y.lim[1], y.lim[2]*1.05, border.col = NA, fill = region.cols[i])
		}
	points(snp.pos, all.nlp, type = "p", pch = 16)
	mtext("-log p value", side = 2, line = 2.5)
	abline(h = -log10(sig.p))
	axis(2)
		
		
	par(mar = c(2, 4, 0,2))
	plot.new()
	all.scores <- as.numeric(all.gene.info[,value])
	y.lim <- c(min(all.scores, na.rm = TRUE), max(all.scores, na.rm = TRUE)*1.05)
	plot.window(xlim = x.lim, ylim = y.lim)
	for(i in 1:nrow(genetic.regions)){
		draw.rectangle(genetic.regions[i,1], genetic.regions[i,2], y.lim[1], y.lim[2]*1.05, border.col = NA, fill = region.cols[i])
		}		
	text(x = as.numeric(all.gene.info[,"gene.position"]), y = as.numeric(all.gene.info[,value]), labels = rownames(all.gene.info), cex = 0.7)
	axis(2)
	mtext(label, side = 2, line = 2.5)	
	if(!is.null(max.fp)){
		abline(h = log10(max.fp)*-1)	
	}

	
	
	
}