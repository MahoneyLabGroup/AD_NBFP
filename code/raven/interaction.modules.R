#This function matches up interactions with functional modules


interaction.modules <- function(fntm.net, Lij, locus.genes, comm.genes){
	
	lgenes.in.net <- lapply(locus.genes, function(x) intersect(V(fntm.net)$entrez.name, x))
	lgenes.names <- lapply(lgenes.in.net, function(x) V(fntm.net)$name[match(x, V(fntm.net)$entrez.name)])

	ints <- which(Lij != 0, arr.ind = TRUE)
	ints[,1] <- rownames(Lij)[as.numeric(ints[,1])]
	ints[,2] <- rownames(Lij)[as.numeric(ints[,2])]	

	all.mod.ints <- vector(mode = "list", length = nrow(ints))

	#====================================================
	#internal functions
	#====================================================	
	direct.link <- function(edges, gene1, gene2){
		edge1.locale <- intersect(which(edges[,1] == gene1), which(edges[,2] == gene2))
		edge2.locale <- intersect(which(edges[,2] == gene1), which(edges[,1] == gene2))
		has.edge <- length(c(edge1.locale, edge2.locale) > 0)
		return(has.edge)
		}
	
	prop.trans <- function(int.mat){
		if(length(int.mat) > 1){
			same.mod <- apply(int.mat, 1, function(x) x[1] == x[2])
			prop.trans <- length(which(!same.mod))/length(same.mod)
			return(prop.trans)		
			}else{
			return("no gene interactions")
			}
		}
	
	num.mod.ints <- function(int.mat){
		if(length(int.mat) > 1){
			u_ints <- unique(int.mat)
			int.counts <- apply(u_ints, 1, function(x) length(intersect(which(int.mat[,1] == x[1]), which(int.mat[,2] == x[2]))))
			final.table <- cbind(u_ints, int.counts)
			colnames(final.table) <- c("mod1", "mod2", "number")
			final.table <- final.table[order(as.numeric(final.table[,3]), decreasing = TRUE),]
			return(final.table)
			}else{
			return("no gene interactions")
			}
		}
	
	#====================================================	
	net.edges <- as_edgelist(fntm.net)
		
	for(i in 1:nrow(ints)){
		locus1.genes <- lgenes.names[[ints[[i,1]]]]
		locus2.genes <- lgenes.names[[ints[[i,2]]]]

		if(length(locus1.genes) > 0 && length(locus2.genes) > 0){ 
			#find all direct edges between these genes
			all.locus.pairs <- cbind(rep(locus1.genes, length(locus2.genes)), rep(locus2.genes, each = length(locus1.genes)))
			connected <- apply(all.locus.pairs, 1, function(x) direct.link(net.edges, x[1], x[2]))
			locus.pairs <- all.locus.pairs[which(connected == 1),]
							
			locus1.modules <- names(comm.genes)[unlist(lapply(locus.pairs[,1], function(x) grep(x, comm.genes)))]
			locus2.modules <- names(comm.genes)[unlist(lapply(locus.pairs[,2], function(x) grep(x, comm.genes)))]
			
			all.mod.ints[[i]] <- cbind(locus1.modules, locus2.modules)
			}else{
			all.mod.ints[[i]] <- "no gene interactions"
			}
		}
	
		all.prop.trans <- lapply(all.mod.ints, prop.trans)
		all.mod.counts <- lapply(all.mod.ints, num.mod.ints)
	
		names(all.prop.trans) <- names(all.mod.counts) <- apply(ints, 1, function(x) paste(x[1], x[2], sep = "-"))
	
		final.result <- list("proportion.trans.module.interactions" = all.prop.trans, "module.interaction.counts" = all.mod.counts)
		return(final.result)
}