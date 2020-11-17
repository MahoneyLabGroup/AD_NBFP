#Given a set of xy coordinates, this function draws a pareto line


pareto.line <- function(x,y, top.percent = 90){
	
	x.top <- get.percentile(x, top.percent)
	y.top <- get.percentile(y, top.percent)
	
	abovex <- which(x >= x.top)
	abovey <- which(y >= y.top)
	
	all.above <- sort(unique(c(abovex, abovey)))
	
	point.coord <- cbind(x[all.above], y[all.above])
	# point.coord <- point.coord[order(point.coord[,1]),]
	
	smoothed <- loess.smooth(point.coord[,1], point.coord[,2])
	
	# plot(x, y)
	# # abline(v = x.top);abline(h = y.top)	
	# # points(point.coord[,1], point.coord[,2], type = "l")
	# points(smoothed$x, smoothed$y, type = "l", col = "red")

	return(smoothed)
	
}