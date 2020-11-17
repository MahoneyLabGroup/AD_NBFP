sort.independently <-
function(mat, margin = 1){
		
	final.mat <- matrix(NA, nrow = dim(mat)[1], ncol = dim(mat)[2])

	if(margin == 1){
		for(i in 1:length(mat[,1])){
			final.mat[i,] <- sort(mat[i,])
			}
		}
	
	if(margin == 2){
		for(i in 1:length(mat[1,])){
			final.mat[,i] <- sort(mat[,i])
		}	
	}	

	return(final.mat)
	
}
