#locus.genes and Lij are optional for highlighting edges between
#nodes that correspond to the locus network
#if color.by.locus is TRUE, this overrides highlight.nodes
#min.edge.weight = 0.15; highlight.nodes = NULL; layout.call = layout_nicely; layout.matrix = NULL; stretch.layout = 1.2; layout.x.adjust = 0.2; layout.y.adjust = 0.1; vertex.size = 10; label.size = 1; locus.genes = NULL; Lij = NULL; show.Lij.only = FALSE; color.by.locus = FALSE; vertex.label = NULL

plot.fntm.net <- function(fntm.net, min.edge.weight = 0.15, highlight.nodes = NULL, layout.call = layout_nicely, layout.matrix = NULL, stretch.layout = 1.2, layout.x.adjust = 0.2, layout.y.adjust = 0.1, vertex.size = 10, label.size = 1, locus.genes = NULL, Lij = NULL, show.Lij.only = FALSE, color.by.locus = FALSE, vertex.label = NULL, heatmap.factor = 1){
	require(pheatmap)
	require(igraph)
	if(color.by.locus && is.null(locus.genes)){stop("I cannot color by locus without locus.genes")}

	low.edges <- which(E(fntm.net)$weight < min.edge.weight)
	if(length(low.edges) > 0){
		trimmed.net <- delete.edges(fntm.net, low.edges)
		}else{
		trimmed.net <- fntm.net
		}

	#delete anything with a 0 degree unless is one
	#of the highlighted nodes
	delete.verts <- which(degree(trimmed.net) == 0)
	not.highlighted <- setdiff(delete.verts, highlight.nodes)
	if(length(not.highlighted) > 0){
		trimmed.net <- delete.vertices(trimmed.net, not.highlighted)
		}
	
	#plot(trimmed.net)
	
	#=====================================================================
	# internal functions
	#=====================================================================	
	
	#This function highlights edges between genes in interacting loci
	highlight.edges <- function(trimmed.net, lgenes.in.net, lgenes.names, all.edges){
		for(i in 1:nrow(all.edges)){
			locus1.genes <- lgenes.names[[all.edges[[i,1]]]]
			locus2.genes <- lgenes.names[[all.edges[[i,2]]]]
			if(length(locus1.genes) > 0 && length(locus2.genes) > 0){
				locus.pairs <- cbind(rep(locus1.genes, length(locus2.genes)), rep(locus2.genes, each = length(locus1.genes)))
				links <- apply(locus.pairs, 1, function(x) shortest_paths(trimmed.net, from = x[1], to = x[2], output = "epath"))
				all.highlight.links <- unique(unlist(lapply(links, function(x) as.vector(x$epath[[1]]))))
				edge.width[all.highlight.links] <- 5
				}
			} #end looping though edges
		return(edge.width)
		} #end highlight.edges()

	#This function trims the network to only the genes and edges
	#in the cape and gene networks
	trim.to.Lij <- function(trimmed.net, lgenes.in.net, lgenes.names, all.edges){
		net.adj <- as.matrix(as_adjacency_matrix(trimmed.net, attr = "weight"))
		new.adj <- matrix(0, nrow = nrow(net.adj), ncol = ncol(net.adj))
		rownames(new.adj) <- colnames(new.adj) <- colnames(net.adj)
		for(i in 1:nrow(all.edges)){
			locus1.genes <- lgenes.names[[all.edges[[i,1]]]]
			locus2.genes <- lgenes.names[[all.edges[[i,2]]]]
			new.adj[locus1.genes, locus2.genes] <- net.adj[locus1.genes, locus2.genes]
			}
		new.net <- graph_from_adjacency_matrix(new.adj, weighted = TRUE)
		V(new.net)$entrez.name <- V(trimmed.net)$entrez.name
		delete.verts <- which(degree(new.net) == 0)
		not.highlighted <- setdiff(delete.verts, highlight.nodes)
		if(length(not.highlighted) > 0){
			new.net <- delete.vertices(new.net, not.highlighted)
			}
		#plot(new.net)		
		return(new.net)
		}
	#=====================================================================	


	#if we want to highlight the edges bewteen loci
	#find all vertices that are in the loci
	#and highlight the edges between genes that match
	#up with the locus network
	edge.width = rep(1, ecount(trimmed.net))
	if(!is.null(locus.genes)){
		if(is.null(Lij)){stop("If locus.genes is set, I also need Lij.")}
		
		edge.width = rep(1, ecount(trimmed.net))
		lgenes.in.net <- lapply(locus.genes, function(x) intersect(V(trimmed.net)$entrez.name, x))
		lgenes.names <- lapply(lgenes.in.net, function(x) V(trimmed.net)$name[match(x, V(trimmed.net)$entrez.name)])
		all.edges <- which(Lij != 0, arr.ind = TRUE)
		all.edges[,1] <- rownames(Lij)[as.numeric(all.edges[,1])]
		all.edges[,2] <- rownames(Lij)[as.numeric(all.edges[,2])]				
		edge.width <- highlight.edges(trimmed.net, lgenes.in.net, lgenes.names, all.edges)
	
	if(show.Lij.only){
		trimmed.net <- trim.to.Lij(trimmed.net, lgenes.in.net, lgenes.names, all.edges)
		# plot(trimmed.net)
		}

		} #end case for highlighting edges


	node.col <- rep("white", vcount(trimmed.net))
	
		if(!is.null(highlight.nodes) && !(color.by.locus)){
			node.locale <- c(match(highlight.nodes, V(trimmed.net)$name), match(as.character(highlight.nodes), V(trimmed.net)$entrez.name))
			node.locale <- node.locale[which(!is.na(node.locale))]
			if(length(node.locale) > 0){
				node.col[node.locale] <- "red"
				}
			}

		if(color.by.locus){
			cols <- categorical_pal(8)
			lgenes.in.net <- lapply(locus.genes, function(x) intersect(V(trimmed.net)$entrez.name, x))
			lgenes.names <- lapply(lgenes.in.net, function(x) V(trimmed.net)$name[match(x, V(trimmed.net)$entrez.name)])
			for(i in 1:length(lgenes.names)){
				node.col[match(lgenes.names[[i]], V(trimmed.net)$name)] <- cols[i]
				}
			}
		
	mypal <- colorRampPalette(c("blue", "#007FFF", "cyan","#7FFF7F", "#7FFF7F", "yellow", "#FF7F00", "red"))
	ColorRamp <- mypal(256)
	possible.weights <- segment.region(0,1,256)
	edge.weights <- E(trimmed.net)$weight
	edge.col <- unlist(lapply(edge.weights, function(x) ColorRamp[get.nearest.pt(possible.weights, x)]))
	
	if(!is.null(layout.matrix)){
		the.layout <- layout.matrix
		}else{
	if(is.null(layout.call)){
		the.layout <- layout_nicely(trimmed.net)
		}else{
		layout.fun <- match.fun(layout.call)	
		the.layout <- layout.fun(trimmed.net)
		}
	}
		
	norm.layout <- norm_coords(the.layout, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
	stretched.layout <- norm.layout*stretch.layout
	stretched.layout[,1] <- stretched.layout[,1] + layout.x.adjust #move along x axis
	stretched.layout[,2] <- stretched.layout[,2] + layout.y.adjust #move along y axis
	
	layout(matrix(c(1,2), nrow = 1), widths = c(1, 0.3))
	if(color.by.locus){
		plot(trimmed.net, layout = stretched.layout, vertex.color = node.col, edge.color = edge.col, vertex.size = vertex.size, rescale = FALSE, label.size = label.size, edge.width = edge.width, vertex.label = vertex.label)		
		legend("topleft", legend = names(locus.genes), fill = cols)
		}else{
		plot(trimmed.net, layout = stretched.layout, vertex.color = node.col, edge.color = edge.col, vertex.size = vertex.size, rescale = FALSE, label.size = label.size, edge.width = edge.width, vertex.label = vertex.label)
		}
	
	barplot(matrix(rep(1, 256), ncol = 1), col = ColorRamp, border = FALSE, axes = FALSE, main = "Edge Weights")
	abline(h = 256*min.edge.weight)
	axis(2, at = segment.region(1,256,6, "ends"), labels = segment.region(0,1,6, "ends"), cex.axis = 2)
	
	layout(matrix(1))
	adj.mat <- as.matrix(as_adjacency_matrix(trimmed.net, attr = "weight"))
	pheatmap(adj.mat^heatmap.factor)
	
	#also sort the adjacency matrix by locus if we have that information
	if(!is.null(locus.genes)){
		lgenes.in.net <- lapply(locus.genes, function(x) intersect(V(trimmed.net)$entrez.name, x))
		lgenes.names <- lapply(lgenes.in.net, function(x) V(trimmed.net)$name[match(x, V(trimmed.net)$entrez.name)])
		lgene.vector <- unlist(lgenes.names)
		lgene.order <- match(lgene.vector, rownames(adj.mat))
		imageWithText(adj.mat[lgene.order,lgene.order], show.text = FALSE, row.names = rownames(adj.mat)[lgene.order], col.names = colnames(adj.mat)[lgene.order])
		idx <- 1
		xv <- c(-1,0)
		yv <- length(lgene.order):0 + 0.5
		par(xpd = TRUE)
		for(i in 1:length(lgenes.names)){
			draw.rectangle(xv[1], xv[2], yv[idx], yv[(idx+length(lgenes.names[[i]]))], fill = categorical_pal(8)[i])
			draw.rectangle(rev(yv)[idx], rev(yv)[(idx+length(lgenes.names[[i]]))], xv[1], xv[2], fill = categorical_pal(8)[i])
			idx <- idx + length(lgenes.names[[i]])
			}
		par(xpd = FALSE)
	}
	
	invisible(trimmed.net)
	
	}