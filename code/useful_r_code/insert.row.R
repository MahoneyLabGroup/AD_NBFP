insert.row <-
function(mat, row, col.or.row = c("col", "row"), after.before = c("after", "before"), location){
		
	
	if(col.or.row == "col"){
		if(after.before == "after"){
			
			#if the new column is to be inserted after the final column
			if(location >= dim(mat)[2]){
				new.mat <- cbind(mat, row) #just paste it on
				}else{ #otherwise, split the matrix and insert the column
					matrix.part1 <- mat[,1:location] #separate out the first chunk of the matrix, all the columns up to the location of insertion
					matrix.part2 <- mat[,(location+1):length(mat[1,])] #separate out the second chunk, all the columns after the location of insertion
					new.mat <- cbind(matrix.part1, row, matrix.part2)
					}
				}
		
		if(after.before == "before"){
			
			#if the new column is to be inserted before the first column
			if(location == 1){
				new.mat <- cbind(row, mat) #just paste it on
				}else{ #otherwise, split the matrix and insert the column			
					matrix.part1 <- mat[,1:(location-1)] #separate out the first chunk of the matrix, all the columns before the location of insertion
					matrix.part2 <- mat[,(location):length(mat[1,])] #separate out the second chunk, all the columns from the location of insertion
					new.mat <- cbind(matrix.part1, row, matrix.part2)	
					}
				}
			}
		
	if(col.or.row == "row"){
		if(after.before == "after"){
			#if the new row is to be inserted after the final row
			if(location >= dim(mat)[2]){
				new.mat <- rbind(mat, row) #just paste it on
				}else{ #otherwise, split the matrix and insert the column
					matrix.part1 <- mat[1:location,] #separate out the first chunk of the matrix, all the columns up to the location of insertion
					matrix.part2 <- mat[(location+1):length(mat[1,]),] #separate out the second chunk, all the columns after the location of insertion
					new.mat <- rbind(matrix.part1, row, matrix.part2)
					}
				}
			
		if(after.before == "before"){
			#if the new row is to be inserted before the first row
			if(location == 1){
				new.mat <- rbind(row, mat) #just paste it on
				}else{ #otherwise, split the matrix and insert the column
					matrix.part1 <- mat[1:(location-1),] #separate out the first chunk of the matrix, all the columns before the location of insertion
					matrix.part2 <- mat[(location):length(mat[1,]),] #separate out the second chunk, all the columns from the location of insertion
					new.mat <- rbind(matrix.part1, row, matrix.part2)	
					}
				}
		}
	
	
	return(new.mat)
	
}
