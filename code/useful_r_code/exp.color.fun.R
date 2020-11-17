exp.color.fun <-
function(x.min, x.max, steepness = 1, num.cols = 256){
	
	x <- seq(0,1,length.out = num.cols)
	y <- exp(steepness*x)-1
	
	#map the 0-1 interval onto the interval from x.min to x.max
	x.range <- x.max - x.min
	scaled.y <- (y*x.range)/max(y)
	

	shifted.y <- scaled.y + x.min
	
	return(shifted.y)
	
	
}
