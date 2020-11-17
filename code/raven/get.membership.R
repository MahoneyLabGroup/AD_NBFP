#This function takes in a network and a tiered matrix
#from iter.cluster2
#and returns a vector of memberships in the same
#order as the network


get.membership <- function(net, tier.matrix, tier.num){
	
	vertex.names <- V(net)$name
    cluster.mem <- tier.matrix[,tier.num]
	mem <- lapply(names(cluster.mem), function(x) strsplit(x, "-")[[1]])
	memV  <- rep(NA, vcount(net))

    for(i in 1:length(mem)){
        gene.locale <- which(vertex.names %in% mem[[i]])
        memV[gene.locale] <- as.numeric(cluster.mem[i])
    }
    return(memV)
}