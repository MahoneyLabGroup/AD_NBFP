#This function pulls the gene network out of a TRiAGE data object

get.net <- function(data.obj){
	
	loci <- names(data.obj$locus.genes)
	genes <- data.obj$cur.gene[loci]
	return(genes)

}