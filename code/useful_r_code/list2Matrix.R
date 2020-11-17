list2Matrix <-
function(listX, label = NULL, preserve.dim = FALSE){
	
	num.elements <- length(listX)

	if(preserve.dim == TRUE){
		#if we want to preserve the dimensions and we do
		#not want a particular piece of the list.
		if(is.null(label)){
		the.matrix <- NULL
		mat.names <- NULL
		for(i in 1:num.elements){
			the.matrix <- rbind(the.matrix, listX[[i]])
			mat.names <- c(mat.names, paste(1:length(listX[[i]][,1]), names(listX)[[i]], sep = ":"))
			}
		rownames(the.matrix) <- mat.names
		}else{
			#if we specified a label to pull out
			list.labels <- names(listX[[1]])
			label.locale <- which(list.labels == label)
			
			if(length(label.locale) == 0){
				stop("The name of your label does not match a label in the list")
				}

			the.matrix <- NULL
			mat.names <- NULL
			for(i in 1:num.elements){
				the.matrix <- rbind(the.matrix, listX[[i]][label.locale][[1]])
				mat.names <- c(mat.names, paste(1:length(listX[[i]][label.locale][[1]][,1]), names(listX)[[i]], sep = ":"))
				}
			rownames(the.matrix) <- mat.names
			}
		}else{


	#use one of these loops if the list is ragged
	if(is.null(label)){
		num.cells <- as.vector(sapply(listX, length))		
		the.matrix <- matrix(NA, nrow = num.elements, ncol = max(num.cells))
		rownames(the.matrix) <- names(listX)
	
		for(i in 1:num.elements){
			if(num.cells[i] > 0){
				the.matrix[i,1:num.cells[i]] <- listX[[i]]
				}
			}
		}else{
			list.labels <- names(listX[[1]])
			label.locale <- which(list.labels == label)

			if(length(label.locale) == 0){
				stop("The name of your label does not match a label in the list")
				}
			num.cells <- as.vector(sapply(listX, function(listX) length(listX[[label.locale]])))

			the.matrix <- matrix(NA, nrow = num.elements, ncol = max(num.cells))
			rownames(the.matrix) <- names(listX)
		
			for(i in 1:num.elements){
				if(num.cells[i] > 0){
					the.matrix[i,1:num.cells[i]] <- listX[[i]][label.locale][[1]]
					}
				}	
			}
		}
		
	return(the.matrix)
	}
