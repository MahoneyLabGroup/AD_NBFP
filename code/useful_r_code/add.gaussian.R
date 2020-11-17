add.gaussian <-
function(v1, percent.total = 5){
	
	gauss.sd <- (percent.total/100)*(1/3)
	noise <- rnorm(length(v1), sd = gauss.sd)
	noisy.vector <- v1+noise
	return(noisy.vector)	
	
}
