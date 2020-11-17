merge.mats <- function(mat.list, keep.all.rows = TRUE){
	
	if(keep.all.rows){
		all.ind <- unique(unlist(lapply(mat.list, rownames)))
		}else{
		all.ind <- Reduce("intersect", lapply(mat.list, rownames))	
		}
		
		all.col <- unlist(lapply(mat.list, colnames))
		
		new.mat <- matrix(NA, nrow = length(all.ind), ncol = length(all.col))
		rownames(new.mat) <- all.ind
		colnames(new.mat) <- all.col
		
		place.col <- function(mat){
			for(i in 1:ncol(mat)){
				from.mat.locale <- match(rownames(new.mat), rownames(mat))
				new.mat[, colnames(mat)[i]] <- mat[from.mat.locale,i]
				}
			return(new.mat)
			}
		
		
		for(i in 1:length(mat.list)){
			new.mat <- place.col(mat.list[[i]])
			}
		
		return(new.mat)

}