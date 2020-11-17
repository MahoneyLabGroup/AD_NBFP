#This function returns nodes of a network that are 
#clustered at the same level after iterative clustering
#by iter.cluster2()
#the module position matrix is returned from iter.cluster2()

clustered.nodes <- function(module.position.matrix, tier = 1){
	
	mod.names <- apply(module.position.matrix[,1:tier,drop=FALSE], 1, 
	function(x) paste(x, collapse = "-"))
	
	u_mod.names <- unique(mod.names)
	
	get.mod.nodes <- function(mod){
		mod.locale <- which(mod.names == mod)
		all.nodes <- rownames(module.position.matrix)[mod.locale]
		split.nodes <- unlist(lapply(all.nodes, function(x) strsplit(x, "-")))
		return(split.nodes)
		}
	
	mod.nodes <- lapply(u_mod.names, get.mod.nodes)
	names(mod.nodes) <- u_mod.names
	return(mod.nodes)
	
	}