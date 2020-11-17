#this function takes in specific genes and returns the subnetwork
#using only those genes


genes2subnet <- function(gene.names, snp.genes){
	

	gene.locale <- which(snp.genes[,"external_gene_name"] %in% gene.names)
	gene.blocks <- snp.genes[gene.locale,"blockID"]
	gene.main <- just.main[gene.blocks,,drop=FALSE]
	
	gene.int <- as.matrix(just.int[gene.blocks, gene.blocks])
	colnames(gene.int) <-  rownames(gene.int) <- snp.genes[gene.locale,"external_gene_name"]
	rownames(gene.main) <- snp.genes[gene.locale,"external_gene_name"]


	#condense the subnetwork down to a unique set of genes
	u_genes <- unique(colnames(gene.int))
	for(g in 1:length(u_genes)){
		gene.idx <- which(colnames(gene.int) == u_genes[g])
		if(length(gene.idx) > 1){
			gene.as.source <- gene.int[gene.idx,,drop=FALSE]
			gene.one.source <- matrix(apply(gene.as.source, 2, function(x) x[which.max(abs(x))]), nrow = 1)
			rownames(gene.one.source) <- u_genes[g]
			colnames(gene.one.source) <- colnames(gene.as.source)
			final.source <- t(unique(cbind(colnames(gene.one.source), t(gene.one.source))))[2,,drop=FALSE]
	
			gene.as.target <- gene.int[, gene.idx,drop=FALSE]
			gene.one.target <- matrix(apply(gene.as.target, 1, function(x) x[which.max(abs(x))]), ncol = 1)
			colnames(gene.one.target) <- u_genes[g]
			rownames(gene.one.target) <- rownames(gene.as.target)
			final.target <- unique(cbind(rownames(gene.one.target), gene.one.target))[,2,drop=FALSE]	

			no.targ.source <- gene.int[-gene.idx, -gene.idx,drop=FALSE]
			add.source <- rbind(final.source[,rownames(no.targ.source),drop=FALSE], no.targ.source)
			add.target <- cbind(final.target[rownames(add.source),,drop=FALSE], add.source)
			target.rows <- rownames(add.target)
			num.target <- apply(add.target, 2, as.numeric)
			rownames(num.target) <- target.rows
		
			#also adjust the main effects
			gene.as.main <- gene.main[gene.idx,,drop=FALSE]
			max.main <- matrix(apply(gene.as.main, 2, function(x) x[which.max(abs(x))]), nrow = 1)
			rownames(max.main) <- u_genes[g]
			colnames(max.main) <- colnames(gene.main)
			pared.main <- rbind(max.main, gene.main[-gene.idx,,drop=FALSE])
			
			gene.int <- num.target
			gene.main <- pared.main
			}
		}#end looping through unique gene names
	
	net.obj <- make.net.obj(int.mat = gene.int, main.mat = gene.main)
	# plot(net.obj, layout = layout_nicely)
	return(net.obj)
	

	
}