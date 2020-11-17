#This function merges QTL from a QTL list 
#such that QTL that overlap are turned into
#a single QTL

merge.qtl <- function(qtl.table){
	
	all.qtl <- unique(c(qtl.table[,1], qtl.table[,2]))
	split.qtl <- strsplit(all.qtl, ":")
		
	#given two positions, this function determines whether they overlap
	does.overlap <- function(qtl1, qtl2){
			
		split.qtl1 <- strsplit(qtl1, ":")
		split.qtl2 <- strsplit(qtl2, ":")
		chr1 <- split.qtl1[[1]][1]
		chr2 <- split.qtl2[[1]][1]		
		if(chr1 != chr2){
			return(FALSE)
			}
		start1 <- as.numeric(split.qtl1[[1]][2]);stop1 <- as.numeric(split.qtl1[[1]][3])
		start2 <- as.numeric(split.qtl2[[1]][2]);stop2 <- as.numeric(split.qtl2[[1]][3])
		region1 <- start1:stop1
		region2 <- start2:stop2
		overlap <- intersect(region1, region2)
		if(length(overlap) > 0){
			return(TRUE)
			}else{
			return(FALSE)
			}	
		}
	
	combine.qtl <- function(qtl.list){
		qtl.chr <- strsplit(qtl.list[1], ":")[[1]][1]
		min.qtl.start <- min(as.numeric(lapply(strsplit(qtl.list, ":"), function(x) x[2])))
		max.qtl.stop <- max(as.numeric(lapply(strsplit(qtl.list, ":"), function(x) x[3])))
		merged.qtl <- paste(qtl.chr, min.qtl.start, max.qtl.stop, sep = ":")
		return(merged.qtl)
		}
	
	cat("Calculating overlaps...\n")
	overlapping.qtl <- matrix(NA, nrow = length(all.qtl), ncol = length(all.qtl))
	rownames(overlapping.qtl) <- colnames(overlapping.qtl) <- all.qtl
	for(i in 1:length(all.qtl)){
		report.progress(i, length(all.qtl))
		overlapping.qtl[i,] <- unlist(lapply(all.qtl, function(x) does.overlap(x, all.qtl[i])))
		}
	diag(overlapping.qtl) <- FALSE
	# image(overlapping.qtl)
	overlapping.qtl[which(overlapping.qtl)] <- 1
	overlapping.qtl[which(!overlapping.qtl)] <- 0
	qtl.net <- graph_from_adjacency_matrix(overlapping.qtl, mode = "undirected")
	# plot(qtl.net)
	comm <- cluster_fast_greedy(qtl.net)$membership
	u_comm <- unique(comm)
	
	new.qtl  <- all.qtl
	for(i in 1:length(u_comm)){
		comm.locale <- which(comm == u_comm[i])
		if(length(comm.locale) > 1){
			new.qtl[comm.locale] <- combine.qtl(all.qtl[comm.locale])
			}
		}
	

	#now replace the qtls in the qtl table with the merged versions
	cat("\nreplacing overlapping QTL with merged QTL...\n")
	for(i in 1:length(new.qtl)){
		all.overlap1 <- unlist(lapply(qtl.table[,1], function(x) does.overlap(new.qtl[i], x)))
		if(length(which(all.overlap1) > 0)){
			qtl.table[which(all.overlap1), 1] <- new.qtl[i]
			}

		all.overlap2 <- unlist(lapply(qtl.table[,2], function(x) does.overlap(new.qtl[i], x)))
		if(length(which(all.overlap2) > 0)){
			qtl.table[which(all.overlap2), 2] <- new.qtl[i]
			}
		}

	return(qtl.table)
	
}