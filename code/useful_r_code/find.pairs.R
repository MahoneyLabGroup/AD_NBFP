#This function takes in two two-column matrices
#It finds the indices in the second matrix
#of the pairs in the first matrix.

find.pairs <- function(mat1, mat2){	
	pair1.text <- apply(mat1, 1, function(x) paste(x[1], x[2], sep = ":"))
	pair2.text <- apply(mat2, 1, function(x) paste(x[1], x[2], sep = ":"))
	pair1.idx <- match(pair1.text, pair2.text)
	# cbind(pair1.text, pair2.text[pair1.idx])
	return(pair1.idx)
	}