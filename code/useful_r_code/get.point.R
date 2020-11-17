#This function gets a point on a line given the slope and intercept


get.point <- function(intercept, slope, x = NULL, y = NULL){
	
	if(is.null(y)){
		point = slope*x+intercept
		}
	if(is.null(x)){
		point = (y-intercept)/slope
		}
	
	return(point)
	
	}