#This function calculates the fitness value of Gij
#it applies a function to the edges dictated by 
#the cape adjacency matrix Jij, for example, sum,
#mean, etc.

fitness <- function(data.obj){
	
	Gij = data.obj$Gij
	Lij = data.obj$Lij
	fun = data.obj$fit.fun
	
	just.int <- as.matrix(Lij[,1:nrow(Lij)])
	
	
	fun <- match.fun(fun)
	fit <- fun(Gij[which(just.int != 0)])
	
	return(fit*-1)
	
	}