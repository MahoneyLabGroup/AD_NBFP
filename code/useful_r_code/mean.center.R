mean.center <-
function(v){
	
	mean.v <- mean(v, na.rm = TRUE)
	centered <- v - mean.v
	return(centered)
	
}
