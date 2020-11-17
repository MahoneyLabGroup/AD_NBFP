#This function builds the locus by locus network based on the genes 
#overlapping the cape network and the svm results
#This uses a qtl table instead of a data object
#qtl.table is the table that comes out of get.genes.in.QTL
#gene.table is an information table for all the genes
#that overlap between the cape network and the SVM.


build_ll_list <- function(qtl.table, gene.list, filter.type = c("entrezgene", "ensembl_gene_id"), mart){
	
	require(biomaRt)
	filter.type <- filter.type[1]
	
	u_qtl <- unique(c(qtl.table[,1], qtl.table[,2]))
	split.qtl <- strsplit(u_qtl, ":")
	qtl.chr <- unlist(lapply(split.qtl, function(x) x[1]))
	qtl.start <- unlist(lapply(split.qtl, function(x) as.numeric(x[2])))
	qtl.end <- unlist(lapply(split.qtl, function(x) as.numeric(x[3])))
		
	#make a gene information table from the overlap genes
	atts <- c("entrezgene", "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")
	overlap.gene.table <- getBM(atts, filter.type, values = gene.list, mart = mart)
	
	genes.by.block <- vector(mode = "list", length = length(u_qtl))
	names(genes.by.block) <- u_qtl
	#for each block get all the genes inside
	for(l in 1:length(qtl.start)){
		chr.locale <- which(overlap.gene.table[,"chromosome_name"] == qtl.chr[l])
		after.start <- which(overlap.gene.table[,"start_position"] >= qtl.start[l])
		before.end <- which(overlap.gene.table[,"end_position"] <= qtl.end[l])
		gene.locale <- Reduce("intersect", list(chr.locale, after.start, before.end))
		genes.by.block[[l]] <- overlap.gene.table[gene.locale,filter.type]
		}
		
	return(genes.by.block)
	
	
}