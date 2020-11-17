#This function plots the adjacency matrix for only 
#true positive genes in the tissue.mat

plot.tp.adj <- function(tissue.mat){
	
	just.tp <- get.tp(tissue.mat)	
	col.order <- match(rownames(just.tp), colnames(just.tp))
	imageWithText(just.tp[,col.order], show.text = FALSE)
		
}