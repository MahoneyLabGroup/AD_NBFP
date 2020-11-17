#This function pulls a network out of the SA results

SA.final.network <- function(data.obj){
	
	
	cur.gene <- data.obj$optimized.output$cur.gene
	final.net <- cur.gene[which(!is.na(cur.gene))]
	return(final.net)
	
	
}