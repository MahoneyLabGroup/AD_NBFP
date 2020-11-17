#This function returns all unique genes in a QTL
#table returned from get.QTL


get.genes.in.QTL <- function(qtl.table, mart){
	
	require(biomaRt)
	atts <- c("external_gene_name", "entrezgene", "ensembl_gene_id", "chromosome_name", "start_position", "end_position")
	fils <- "chromosomal_region"

	all.qtl <- unique(c(qtl.table[,1], qtl.table[,2]))
	all.genes <- getBM(atts, fils, values = all.qtl, mart = mart)
	return(all.genes)
	
	
}