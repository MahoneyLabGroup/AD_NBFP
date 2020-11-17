z.score <-
function(v){
	
	z <- (v - mean(v))/sd(v)

	return(z)
	
}
