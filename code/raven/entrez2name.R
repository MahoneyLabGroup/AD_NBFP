#This function takes in a matrix of ensembl gene
#ids and translates them to gene names
#it is used in get.gene.stats()

entrez2name <- function(entrez.mat, mart){
	
	u_ids <- unique(as.vector(entrez.mat[which(!is.na(entrez.mat))]))
	id.info <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", u_ids, mart)
	final.mat <- entrez.mat
	
	for(i in 1:length(u_ids)){
		id.locale <- which(id.info[,1] == u_ids[i])[1]
		test <- which(entrez.mat == u_ids[i])
		final.mat[which(entrez.mat == u_ids[i])] <- id.info[id.locale,2]
		}
	
	return(final.mat)
	
}