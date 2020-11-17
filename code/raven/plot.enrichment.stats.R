#plot enrichment stats for modules

plot.enrichment.stats <- function(path = "."){
		
	all.limits <- list.files(pattern = "enrichment.for.modules")	
	all.mod.size <- sort(as.numeric(unlist(lapply(strsplit(all.limits, "\\."), function(x) x[4]))))
	
	# plot individual stats
	#look at the overlap in terms between modules
	mod.pair.stats <- vector(mode = "list", length = length(all.limits))
	for(e in 1:length(all.limits)){
		all.enrich <- readRDS(paste0("enrichment.for.modules.", all.mod.size[e], ".RData"))
		mods <- readRDS(paste0("Module.Membership.", all.mod.size[e], ".RData"))
		u_mods <- unique(mods)
		all.mod.pairs <- matrix(NA, nrow = length(all.enrich), ncol = length(all.enrich))
		rownames(all.mod.pairs) <- colnames(all.mod.pairs) <- u_mods
		for(i in 1:length(all.enrich)){
			for(j in 1:length(all.enrich)){
				all.mod.pairs[i,j] <- jaccard.ind(all.enrich[[i]][,"term.name"], all.enrich[[j]][,"term.name"])	
			}
		}
		all.mod.pairs <- signif(all.mod.pairs, 2)
		mod.pair.stats[[e]] <- all.mod.pairs
		pdf(paste0("Mod.Pairs.", all.mod.size[e], ".pdf"))
		a <- pheatmap(all.mod.pairs)
		imageWithText(all.mod.pairs[a$tree_row$order, a$tree_col$order])
		n.terms <- unlist(lapply(all.enrich, function(x) length(unique(x[,"term.name"]))))
		n.genes <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
		barplot(n.terms[a$tree_col$order], las = 2, main = "Number of terms in each module")
		barplot(n.genes[a$tree_col$order], las = 2, main = "Number of genes in each module")
		plot(n.terms, n.genes)	
		dev.off()
		} #end looping through max mod size


	#plot group stats
	all.n <- vector(mode = "list", length = length(all.limits))
	for(i in 1:length(all.limits)){
		all.enrich <- readRDS(paste0("enrichment.for.modules.", all.mod.size[i], ".RData"))
		mods <- readRDS(paste0("Module.Membership.", all.mod.size[i], ".RData"))
		n.terms <- unlist(lapply(all.enrich, function(x) length(unique(x[,"term.name"]))))
		u_mods <- sort(unique(mods))
		n.genes <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
		all.n[[i]] <- cbind(n.terms, n.genes)
		}
		
	all.x <- unlist(lapply(all.n, function(x) x[,1]))
	all.y <- unlist(lapply(all.n, function(x) x[,2]))

	pdf("Bulk.Module.Stats.pdf")
	plot.new()
	plot.window(xlim = c(min(all.x, na.rm = TRUE), max(all.x, na.rm = TRUE)), ylim = c(min(all.y, na.rm = TRUE), max(all.y, na.rm = TRUE)))	
	for(i in 1:length(all.n)){
		mat.order <- order(as.numeric(rownames(mod.pair.stats[[i]])))
		con.mat <- mod.pair.stats[[i]][mat.order, mat.order]
		for(j in 1:(nrow(mod.pair.stats[[i]])-1)){
			for(k in (j+1):ncol(mod.pair.stats[[i]])){
				segments(all.n[[i]][j,1], all.n[[i]][j,2], all.n[[i]][k,1], all.n[[i]][k,2], lwd = (con.mat[j,k]), col = "gray")
				}
			}
		points(all.n[[i]][,1], all.n[[i]][,2], col = i, pch = 16)

		}
	axis(1); axis(2)
	legend("topleft", col = 1:length(all.mod.size), pch = 16, legend = all.mod.size)

	off.diag <- lapply(mod.pair.stats, function(x) x[upper.tri(x)])
	boxplot(off.diag, names = all.mod.size, main = "Jaccard Indices")
	dev.off()	
	
}