#This function finds the slope between adjacent points
#in x y vectors


instant.slope <- function(x, y){
	
	consec.mat <- consec.pairs(1:length(x))

	slope <- function(x0, x1, y0, y1){
		return((y1-y0)/(x1-x0))
		}
	
	all.slopes <- apply(consec.mat, 1, function(v) slope(x[v[1]], x[v[2]], y[v[1]], y[v[2]]))
	
	return(all.slopes)

}