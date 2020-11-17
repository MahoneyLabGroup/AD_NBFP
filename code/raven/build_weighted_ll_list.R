#This function builds the locus by locus network based on the genes 
#overlapping the cape network and the svm results
#This uses a qtl table instead of a data object
#qtl.table is the table that comes out of get.genes.in.QTL
#gene.table is an information table for all the genes
#that overlap between the cape network and the SVM.
#This function differs from build_ll_list() in that
#it assigns weights to each gene in the loci from a second
#list of trait.related genes. This is a soft filter on the
#list that allows us to consider trait-related genes more
#highly than non-trait-related genes, but doesn't take
#the trait-related genes out of the running completely
#if trait.genes is NULL, all genes get the same weight
#otherwise, the genes in trait.genes will be weighted
#higher than other locus genes based on the weight.ratio

build_weighted_ll_list <- function(qtl.table, gene.list, trait.genes = NULL, weight.ratio = 2, filter.type = c("entrezgene", "ensembl_gene_id"), mart){
	
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
	
	weighted.blocks <- lapply(genes.by.block, function(x) cbind(x, rep(1, length(x))))
	
	#now assign weights based on the trait gene list
	if(!is.null(trait.genes)){
		trait.gene.locale <- lapply(genes.by.block, function(x) which(x %in% trait.genes))
		for(i in 1:length(trait.gene.locale)){
			weighted.blocks[[i]][trait.gene.locale[[i]],2] <- weight.ratio
			weighted.blocks[[i]][,2] <- as.numeric(weighted.blocks[[i]][,2])	/max(as.numeric(weighted.blocks[[i]][,2]))
			}
		}
		
	return(weighted.blocks)
	
}