#This function returns just the true positive matrix

get.tp <- function(tissue.mat){
	
	tp.genes <- row.names(tissue.mat)
	col.locale <- which(dimnames(tissue.mat)[[2]] %in% tp.genes)
	
	just.tp <- as.matrix(tissue.mat[,col.locale])
	return(just.tp)
	}