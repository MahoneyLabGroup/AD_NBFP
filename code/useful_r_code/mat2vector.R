#This function turns a matrix into a vector either by row or by column



	mat2vector <- function(mat, by.row = TRUE){
		
		
		if(by.row){
			new.v <- unlist(apply(mat, 1, list))
			}else{
			new.v <- unlist(apply(mat, 2, list))
			}

		return(new.v)
		
			
		}
