#This function returns the slope between two points

get.slope <- function(x0, x1=x0, y0, y1=y0){
	
	slope <- (y1-y0)/(x1-x0)
	return(slope)
	
	
}