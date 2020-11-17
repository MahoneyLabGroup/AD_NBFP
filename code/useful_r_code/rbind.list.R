rbind.list <-
function(listX){

		mat <- NULL
		col.labels <- NULL
		for(i in 1:length(listX)){
			mat <- rbind(mat, listX[[i]])
			if(length(rownames(listX[[i]])) == 0){
				col.labels <- c(col.labels, names(listX)[i])	
				}else{
					col.labels <- c(col.labels, rownames(listX[[i]]))
					}
			}
			
		mat <- matrix(mat, dim(mat)[1], dim(mat)[2])
		rownames(mat) <- col.labels
		
		return(mat)
		
	}
