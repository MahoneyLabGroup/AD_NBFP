#This function reads in a set of NetWAS files
#and sorts genes based on their SVM score 
#across all traits

Netwas_to_gene_list <- function(traits = c("ET1", "ET2"), path = "."){
	
	cur.dir <- getwd()
	setwd(path)
	
	all.files <- list.files(pattern = "csv")
	trait.locale <- unlist(lapply(traits, function(x) grep(x, all.files)))
	
	all.tables <- lapply(all.files[trait.locale], function(x) read.table(x, sep = "\t", header = TRUE, stringsAsFactors = FALSE))	

	u_genes <- unique(unlist(lapply(all.tables, function(x) x[,1])))

	svm.table <- matrix(NA, nrow = length(u_genes), ncol = length(traits))
	for(i in 1:length(all.tables)){
		gene.locale <- match(u_genes, all.tables[[i]][,1])
		# head(cbind(u_genes, all.tables[[i]][gene.locale,1]))
		svm.table[,i] <- all.tables[[i]][gene.locale,3]
		}
	rownames(svm.table) <- u_genes
	colnames(svm.table) <- traits

	max.svm <- 	apply(svm.table, 1, max)
	ordered.table <- svm.table[order(max.svm, decreasing = TRUE),]

	return(ordered.table)
	
}