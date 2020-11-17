
plot.with.col <- function(x, y, z, low.col = "blue", high.col = "red", breaks = 50, xlab = "x axis", ylab = "y axis", zlab = "z label", main = "title"){
	
	l.col <- get.color(low.col)[2]
	h.col <- get.color(high.col)[2]

	rbPal <- colorRampPalette(c(l.col,h.col))
	dat.cols <- rbPal(breaks)[as.numeric(cut(z,breaks = breaks))]
	seq.int <- (max(z, na.rm = TRUE) - min(z, na.rm = TRUE))/breaks
	col.mat <- seq(min(z, na.rm = TRUE), max(z, na.rm = TRUE), seq.int)
	
	layout.mat <- matrix(c(1,2), ncol = 2)
	layout(layout.mat, widths = c(1,0.2))
	par(mar = c(5,4,4,2))
	plot(x, y, pch = 20,col = dat.cols, xlab = xlab, ylab = ylab)
	par(mar = c(3,2.5,5,2))	
	image(1, col.mat, matrix(data=col.mat, ncol=length(col.mat),nrow=1), col=rbPal(breaks), xlab="",ylab="",xaxt="n", cex.axis = 2, main = zlab)
	
	
}