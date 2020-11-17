plot.example.net <- function(Lij, gene.ids, mart, plot.label = NULL, mark.neighborhood = NULL, node.size.by.degree = TRUE, layout.call = layout_nicely, layout.matrix = NULL, stretch.layout = 1.2, layout.x.adjust = 0.2, layout.y.adjust = 0.1, manual.layout = FALSE, show.node.names = TRUE, vertex.alpha = 1, vertex.col = "lightblue", highlight.col = "purple", positive.edge.col = "brown", negative.edge.col = "blue", edge.curve = 0.1){
	
	require(igraph)

	if(!is.null(gene.ids)){
		gene.names <- getBM(c("entrezgene", "external_gene_name"), "entrezgene", gene.ids, mart)
		names.by.locus <- gene.ids
		for(i in 1:length(names.by.locus)){
			gene.locale <- which(gene.names[,1] == gene.ids[i])
			if(length(gene.locale) > 0){
				names.by.locus[i] <- gene.names[gene.locale,2]	
				}
			}
		na.locale <- which(is.na(names.by.locus))
		names.by.locus[na.locale] <- names(names.by.locus)[na.locale]
		adj.mat <- Lij
		adj.mat[which(is.na(adj.mat))] <- 0
		
		locus.locale <- match(names(names.by.locus), colnames(adj.mat))
		if(length(locus.locale) > 0){
			adj.mat <- adj.mat[locus.locale, locus.locale]
			colnames(adj.mat) <- rownames(adj.mat) <- names.by.locus
			}

		}else{
		adj.mat <- Lij
		adj.mat[which(is.na(adj.mat))] <- 0
		colnames(adj.mat) <- rownames(adj.mat) <- NULL

		}


	net <- graph.adjacency(adj.mat, mode = "directed", weighted = TRUE)
	e.col <- rep(positive.edge.col, ecount(net))
	e.col[which(E(net)$weight < 0)] <- negative.edge.col
	
	E(net)$arrow.size <- 0.4
	E(net)$width <- 1+abs(E(net)$weight)/12

	if(!is.null(layout.matrix)){
		the.layout <- layout.matrix
		}else{
		if(is.null(layout.call)){
			the.layout <- layout_nicely(net)
			}else{
			layout.fun <- match.fun(layout.call)	
			the.layout <- layout.fun(net)
			}
		}
		
	norm.layout <- norm_coords(the.layout, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
	stretched.layout <- norm.layout*stretch.layout
	stretched.layout[,1] <- stretched.layout[,1] + layout.x.adjust #move along x axis
	stretched.layout[,2] <- stretched.layout[,2] + layout.y.adjust #move along y axis

	vcol <- rep(adjustcolor(vertex.col, alpha = vertex.alpha), vcount(net))
	if(length(mark.neighborhood) > 0){
		igf.interactors <- unique(c(names(which(adj.mat[mark.neighborhood,] != 0)), names(which(adj.mat[,mark.neighborhood] != 0))))
		to.mark <- c(mark.neighborhood, igf.interactors)	
		vcol[match(to.mark, V(net)$name)] <- adjustcolor(highlight.col, alpha = vertex.alpha)
		}
	V(net)$color <- vcol
	
	if(node.size.by.degree){
		modes = c("all", "out", "in")
		m = 1
		# for(m in 1:length(modes)){
			deg <- degree(net, mode=modes[m])
			V(net)$size <- deg*3
			}

		# quartz(width = 8, height= 8)
		par(xpd = TRUE)
		if(show.node.names){
		plot(net, edge.curved = edge.curve, layout = stretched.layout, vertex.frame.color = "#ffffff", edge.color = e.col, rescale = FALSE, vertex.shape = "rectangle", main = plot.label)
		}else{
		plot(net, edge.curved = edge.curve, layout = stretched.layout, vertex.frame.color = "#ffffff", edge.color = e.col, rescale = FALSE, vertex.label = NA, vertex.shape = "rectangle", main = plot.label)
			}
		par(xpd = FALSE)
		# }


	
	if(manual.layout){
		tkp.id <- tkplot(net, layout = stretched.layout)
		done <- readline(prompt = "Press return when ready:\n")
		the.layout <- tkplot.getcoords(tkp.id)
		tkplot.close(tkp.id)
		# saveRDS(layout.matrix, "layout.matrix.RData")
		}
		
		invisible(the.layout)

	
}