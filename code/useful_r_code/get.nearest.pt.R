get.nearest.pt <- function(v, pt){
	# plot(abs(pt - v))
	return(which.min(abs(v - pt)))
	}