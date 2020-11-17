#This function clusters a network iteratively with fast greedy
#clustering. It uses module size as a stopping criterion
#to suppress iterative clustering, set max.mod.size to NULL

iter.cluster <- function(net, max.mod.size = 200){

	mods <- cluster_fast_greedy(net)$membership
	
	if(is.null(max.mod.size)){
		return(mods)
		}
	
	u_mods = sort(unique(mods))			
	mod.size <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
			
	cluster.large.mods <- function(mods){
		u_mods = sort(unique(mods))
		mod.size <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
		large.mods <- which(mod.size > max.mod.size)
		for(i in u_mods[large.mods]){
			large.mod.locale <- which(mods == i)
			subnet <- induced_subgraph(net, vids = large.mod.locale)
			sub.mods <- cluster_fast_greedy(subnet)$membership
			u_sub.mods <- sort(unique(sub.mods))
			mods[large.mod.locale] <- sub.mods + max(mods)
			}
		return(mods)
		}
	
	
	while(any(mod.size > max.mod.size)){
		mods <- cluster.large.mods(mods)
		u_mods <- sort(unique(mods))
		mod.size <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
		# print(mod.size)
		}


	#change module numbers such that they are consecutive
	new.mods <- rep(NA, length(mods))
	mod.name <- 1
	mod.num <- 1
	
	for(mod.num in 1:max(mods)){
		mod.locale <- which(mods == mod.num)
		if(length(mod.locale) > 0){
			new.mods[mod.locale] <- mod.name
			mod.name <- mod.name + 1
			}
		}

	return(new.mods)


}





