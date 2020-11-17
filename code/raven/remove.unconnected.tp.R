#This function prunes a tissue adjacency matrix 
#to remove genes with rows and columns that don't
#vary. It also sorts the matrix such that all annotated
#positives are in the first section of the matrix

remove.unconnected.tp <- function(tissue.mat){
	
	tp.genes <- row.names(tissue.mat)
	tn.genes <- setdiff(colnames(tissue.mat), rownames(tissue.mat))
	
	tp.locale <- match(tp.genes, colnames(tissue.mat))
	tn.locale <- match(tn.genes, colnames(tissue.mat))
	
	just.tp <- as.matrix(tissue.mat[,tp.locale])
	just.tn <- as.matrix(tissue.mat[,tn.locale])
			
	row.means <- rowSums(just.tp)
	zero.locale <- which(row.means == 0)
	
	
	if(length(zero.locale) > 0){
		cat("Removing", length(zero.locale), "genes with no connections.\n")
		just.tp <- just.tp[-zero.locale,]			
		just.tp <- just.tp[,-zero.locale]
		just.tn <- just.tn[-zero.locale,]
		}
		
	full.mat <- cbind(just.tp, just.tn)

	return(full.mat)
		
}