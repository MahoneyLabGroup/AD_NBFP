#This function take a matrix, and adds
#the minimum value to the diagonal that 
#will make the matrix positive definite.


pos.def <- function(matX, max.eig = 0.1){	
	
	min.eig <- min(eigen(matX)$values)
	if(min.eig < 0){	
		diag(matX) = diag(matX) - min.eig
		}
		
	return(matX)
	}