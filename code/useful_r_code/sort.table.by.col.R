sort.table.by.col <- function(a, b = NULL, sort.col = 1, decreasing = FALSE){
	
	if(length(dim(a)) > 0){
		#then we have a matrix. Sort by the sort column
		sorted <- sort(a[,sort.col], index.return = TRUE, decreasing = decreasing)
		
		sorted.a <- matrix(NA, ncol = length(a[1,]), nrow = length(a[,1]))

		for(i in 1:length(a[1,])){
			sorted.a[,i] <- a[sorted$ix, i]
			}
		
		}else{
			#if we just put in two vectors, sort the second one relative to the first
			#and return a matrix of the two vectors
			
			sorted.v <- sort(a, index.return = TRUE)
			sorted.b <- b[sorted.v$ix]
			sorted.a <- matrix(c(sorted.v$x, sorted.b), ncol = 2, byrow = FALSE)
			}
	
		
	dimnames(sorted.a) <- dimnames(a)
	return(sorted.a)
	
	}
