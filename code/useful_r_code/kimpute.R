#This function uses k nearest neighbors to impute missing values
#columns with more than missing.thresh missing data are removed
#entirely

kimpute <- function(mat, k = 10, missing.thresh = 10){

		ind.missing <- which(apply(mat, 1, function(x) length(which(is.na(x)))) > 0)
		
		#go through the individuals, find the k nearest neighbors, 
		#and vote on each of their missing data point
		if(length(ind.missing) > 0){
			# heatmap(chr.geno, Rowv = NA)
			data.dist <- as.matrix(dist(mat))
			for(i in 1:length(ind.missing)){
				ind <- mat[ind.missing[i],,drop=FALSE]
				nearest.neighbors <- as.numeric(names(sort(data.dist[ind.missing[i],])[1:k]))
				neighbors <- mat[nearest.neighbors,]
				missing <- which(is.na(ind))
				ind[missing] <- apply(neighbors, 2, function(x) mean(x, na.rm = TRUE))[missing]
				mat[ind.missing[i],] <- ind
				}
			}

		#now take out any individuals with more than the threshold of missing data
		ind.missing <- which(apply(mat, 1, function(x) length(which(is.na(x)))) > 0)
		perc.missing <- apply(mat, 1, function(x) length(which(is.na(x)))/length(x))*100
		lots.missing <- which(perc.missing > missing.thresh)

		if(length(lots.missing) > 0){
			cat("Removing", length(lots.missing), "rows due to missing data.\n")
			mat <- mat[-lots.missing,,drop=FALSE]
			}
	
	return(list("imputed.matrix" = mat, "rows.removed" = lots.missing))
}