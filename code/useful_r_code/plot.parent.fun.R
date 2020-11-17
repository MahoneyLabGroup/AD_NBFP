#This function shows a network describing a set
#of functions and their parent functions
#all functions with foreach loops
#library(cape)
#library(mvbutils)

#foreach.fun <- c("calc.p", "error.prop", "filter.hwe", "genome.wide.threshold.1D.parallel", "get.pairs.for.pairscan", "impute.missing.geno", "kinship", "one.pairscan.parallel", "one.singlescan", "pairscan.kin", "singlescan")
# cape.fun <- find.funs("package:cape")
# sub.fun are the functions for which you want the parents plotted
# all.fun are all the function names in the package
# plot.parent.fun(foreach.fun, cape.fun)
plot.parent.fun <- function(sub.fun, all.fun){
	
	library(mvbutils)
	
	#get the dependency matrix for all functions
	dep.mat <- foodweb(all.fun, plotting = FALSE)$funmat

	# library(igraph)
	# net <- graph.adjacency(dep.mat, mode = "directed")
	# plot(net)
	
	#we need to find all parents to each function with a foreach loop
	get.parents <- function(fun.names){
		all.parents <- NULL
		for(i in 1:length(fun.names)){
			#iterate through the parents to trace back to all parents
			#find the par.fun among the child names
			child.fun.locale <- which(colnames(dep.mat) == fun.names[i])
			#get all parents for par.fun
			parents <- rownames(dep.mat)[which(dep.mat[,child.fun.locale] == 1)]
			all.parents <- c(all.parents, parents)
			}
		return(all.parents)
		}
	
	ancestors  <- list()
	iter <- 1
	parents <- get.parents(sub.fun)
	while(length(parents) > 0){
		ancestors[[iter]] <- unique(parents)
		parents <- get.parents(parents)
		iter = iter + 1
		}
		
	#now look at how everyone is connected
	foodweb(c(unlist(ancestors), sub.fun))
	return(ancestors)
}