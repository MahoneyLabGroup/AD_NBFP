#This function clusters a network iteratively with fast greedy
#clustering. It uses module size as a stopping criterion.
#To suppress iterative clustering, set max.mod.size to NULL
#If a min.mod.size is set, modules smaller than the minimum
#size are re-grouped with the cluster immediately above them
#in the hierarchy. This step can re-generate clusters larger
#than the max, which should be taken into account when setting
#both the min and max.

iter.cluster2 <- function(net, max.mod.size = 400, min.mod.size = 100, plot.cluster.size = FALSE){
	if(is.null(V(net)$name)){V(net)$name <- 1:vcount(net)}

	clust <- cluster_fast_greedy(net)
	mods <- clust$membership
	
	if(is.null(max.mod.size)){
		return(mods)
		}

	#==========================================
	# internal functions
	#==========================================
	make.clust.tree <- function(clust.obj){
		clust.tree <- vector(mode = "list", length = length(clust.obj))
		for(i in 1:length(clust.obj)){
			clust.tree[[i]] <- clust.obj[[i]]
			}
		return(clust.tree)
		}
				
	get.mod.size <- function(clust.branch){
		if(class(clust.branch) != "list"){return(length(clust.branch))}
		if(class(clust.branch) == "list"){lapply(clust.branch, get.mod.size)}	
		}	
	
	subcluster <- function(clust.branch){
		if(class(clust.branch) != "list"){
			mod.size <- length(clust.branch)
			if(mod.size > max.mod.size){
				subnet <- induced_subgraph(net, vids = clust.branch)
				V(subnet)$name <- clust.branch
				sub.clust <- cluster_fast_greedy(subnet)
				sub.branch <- make.clust.tree(sub.clust)
				
				#check the new cluster sizes
				#if any are below the minimum, merge the
				#the smallest branches together until none
				#are below the minimum
				branch.sizes <- sapply(sub.branch, length)
				num.branches <- length(branch.sizes)
				iter = 1
				while(any(branch.sizes < min.mod.size) && iter < num.branches){
					test.clust <- as.hclust(sub.clust)
					cut.clust <- cutree(test.clust, k = (num.branches - iter))
					sub.branch <- vector(mode = "list", length = (num.branches - iter))
					for(i in 1:length(sub.branch)){
						sub.branch[[i]] <- names(cut.clust)[which(cut.clust == i)]
					}
					branch.sizes <- sapply(sub.branch, length)
					iter = iter + 1
				}
				return(sub.branch)
				}else{
				return(clust.branch)	
				}
			}
			if(class(clust.branch) == "list"){
				lapply(clust.branch, subcluster)
				}
			}	
			
		assign.mod.names <- function(clust.branch){
			if(class(clust.branch) != "list"){
				clust.name <- paste(clust.branch, collapse = "-")
				return(clust.name)
				}
			if(class(clust.branch) == "list"){
				lapply(clust.branch, assign.mod.names)
				}
			}	
			
	#==========================================
	
		
	final.clust <- make.clust.tree(clust)
	all.mod.size <- unlist(get.mod.size(final.clust))

	while(any(all.mod.size > max.mod.size)){
		final.clust <- subcluster(final.clust)
		all.mod.size <- unlist(get.mod.size(final.clust))
		}
	
	mod.name.tree <- assign.mod.names(final.clust)
	mod.names <- unlist(mod.name.tree)
	mod.num <- 1:length(mod.names)

	mod.pos <- sapply(mod.names, function(x) grep(x, mod.name.tree))
	clust.name.tree <- unlist(mod.name.tree, recursive = FALSE)
	while(class(clust.name.tree) == "list"){
		mod.pos <- cbind(mod.pos, sapply(mod.names, function(x) grep(x, clust.name.tree)))
		clust.name.tree <- unlist(clust.name.tree, recursive = FALSE)	
		}

	if(class(mod.pos) == "integer"){
		mod.pos.mat <- matrix(mod.pos, ncol = 1)
		rownames(mod.pos.mat) <- names(mod.pos)
		mod.pos <- mod.pos.mat
	}
	colnames(mod.pos) <- paste("tier", 1:ncol(mod.pos), sep = "")

	return(mod.pos)

	}





