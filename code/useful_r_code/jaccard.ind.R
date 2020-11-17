#This function calculates the jaccard index between two vectors

jaccard.ind <- function(v1, v2){
	num <- intersect(v1, v2)
	denom <- union(v1, v2)
	jc <- length(num)/length(denom)
	return(jc)
	}