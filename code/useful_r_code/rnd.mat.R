#This function makes a random symmetric matrix


rnd.mat <- function(n.row, percent.filled){
	
	mat <- matrix(0, n.row, n.row)
	num <- ceiling((length(mat)*(percent.filled/100))/2)
	rnd.num <- runif(num)
	all.coord <- which(upper.tri(mat), arr.ind = TRUE)
	rnd.pos <- sample(1:nrow(all.coord), num)
	for(i in 1:length(rnd.pos)){
		mat[all.coord[rnd.pos[i],1], all.coord[rnd.pos[i],2]] <- rnd.num[i]
		}

	# image(mat)	
	final.mat <- mat + t(mat)
	# quartz();image(final.mat)
	return(final.mat)
	
}
