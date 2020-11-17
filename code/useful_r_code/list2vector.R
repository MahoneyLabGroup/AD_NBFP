list2vector <-
function(listX, list.element = NULL){
	vectorX <- NULL
	num.elements <- length(listX)

	for(i in 1:num.elements){
		if(is.null(list.element)){
			vectorX <- c(vectorX, listX[[i]])
			}else{
				vectorX <- c(vectorX, listX[[i]][list.element])
				}
		}

	return(vectorX)
	
	}
