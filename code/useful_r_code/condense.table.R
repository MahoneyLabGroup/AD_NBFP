#This function condenses a table to unique rows
#and concatenates the distinct information from 
#remaining rows. This works for tables returned
#from biomaRt where distinct information in one
#column generates individual rows, so a lot of 
#information is repeated.


condense.table <- function(tableX, condense.by, col.to.collapse, col.to.concat){
	
	keep.col <- c(condense.by, col.to.collapse, col.to.concat)
	
	u_entries <- unique(tableX[which(!is.na(tableX[,condense.by])),condense.by])
	new.table <- matrix(NA, nrow = length(u_entries), ncol = length(keep.col))
	colnames(new.table) <- colnames(tableX)[keep.col]
	
	for(i in 1:length(u_entries)){
		entry.locale <- which(tableX[,condense.by] == u_entries[i])
		to.cat <- tableX[entry.locale,col.to.concat,drop=FALSE]
		cat.col <- unlist(apply(to.cat, 2, function(x) paste(unique(x), collapse = "; ")))
		new.table[i,] <- unlist(c(u_entries[i], tableX[entry.locale[1],col.to.collapse], cat.col))	
		}
	
	return(new.table)
	
}