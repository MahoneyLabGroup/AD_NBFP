#This function plots a group of enrichments in a heatmap
#of -log p values
#max.char limits the number of characters shown for a term's name
#which helps keep the heatmap visible even when term names are long.
#if max.char is NULL, names are not trimmed.
#transformation is a function applied to the final matrix to help with
#visualization, such as log or sqrt

plot.enrichment.group <- function(enrichment.list, n.terms = 10, max.char = 40, 
cluster_cols = TRUE, cluster_rows = TRUE, transformation = NULL, plot.results = TRUE,
plot.label = NULL, sort.by = c("p.value", "default")){
	
	if(length(enrichment.list) == 1){
		term.mat <- plot.enrichment.vis(enrichment.list[[1]], num.terms = n.terms, order.by = sort.by, 
		plot.label = plot.label)
		}else{

		sort.by <- sort.by[1]

		if(sort.by == "p.value"){
			sorted.enrich <- lapply(enrichment.list, function(x) x[order(x[,"p.value"], decreasing = FALSE),])
		}else{
			sorted.enrich <- enrichment.list
		}
		
		get.terms <- function(enrich, n.terms){
			n.row <- nrow(enrich)
			if(n.row == 0){return(NULL)}
			keep <- min(c(n.row, n.terms))
			term.vals <- enrich[1:keep,c("term.name", "p.value")]
			term.vals[,2] <- -log10(term.vals[,2])
			return(term.vals)
			}
			
		trim.name <- function(term.name, max.char){
			split.name <- strsplit(term.name, "")
			name.len <- length(split.name[[1]])
			if(name.len < max.char){
				return(term.name)
				}else{
				paste.name <- paste0(paste0(split.name[[1]][1:max.char], collapse = ""), "...", collapse = "")
				return(paste.name)
				}	
			}

		all.terms <- lapply(sorted.enrich, function(x) get.terms(x, n.terms))
		u_terms <- unique(unlist(lapply(all.terms, function(x) x[,1])))
		
		term.mat <- matrix(0, nrow = length(u_terms), ncol = length(enrichment.list))
		rownames(term.mat) <- u_terms
		colnames(term.mat) <- names(enrichment.list)
		for(i in 1:length(enrichment.list)){
			term.idx <- match(all.terms[[i]][,1], rownames(term.mat))
			term.mat[term.idx,i] <- all.terms[[i]][,2]
			}
			
		if(!is.null(max.char)){
			trimmed.names <- sapply(rownames(term.mat), function(x) trim.name(x, max.char))
			rownames(term.mat) <- trimmed.names
			}
			
		if(!is.null(transformation)){
			trans.fun <- match.fun(transformation)
			term.mat <- trans.fun(term.mat)	
			}

		if(plot.results){
			if(!is.null(plot.label)){
				pheatmap(term.mat, cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
				main = plot.label)
			}else{
				pheatmap(term.mat, cluster_cols = cluster_cols, cluster_rows = cluster_rows)
			}
		}
		}
		invisible(term.mat)
	}