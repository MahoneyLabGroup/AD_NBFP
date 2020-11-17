#This function adjusts a matrix (matX) of values based on another
#matrix of values (adj.mat). Because residuals() doesn't return values
#for NA's, this function keeps track of NA positions and places
#The corrected values in the right place.

adjust <- function(matX, adj.mat){

	adj.na <- which(is.na(adj.mat), arr.ind = TRUE)
	adj.na.idx <- unique(adj.na[,1])
	
	new.mat <- matrix(NA, nrow = nrow(matX), ncol = ncol(matX))
	for(i in 1:ncol(new.mat)){
		na.locale <- union(which(is.na(matX[,i])), adj.na.idx)
		not.na.locale <- setdiff(1:nrow(matX), na.locale)
		adjV <- residuals(lm(matX[,i]~adj.mat))
		new.mat[not.na.locale,i] <- adjV		
		}

	colnames(new.mat) <- colnames(matX)
	rownames(new.mat) <- rownames(matX)

	return(new.mat)
}
