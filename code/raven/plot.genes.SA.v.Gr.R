plot.genes.SA.v.Gr <- function(rank.table.SA, id.table.SA, rank.table.Gr, id.table.Gr){
	
	
	all.genes <- unique(c(as.vector(id.table.SA), as.vector(id.table.Gr)))
	all.genes <- all.genes[which(!is.na(all.genes))]
	
	num.table <- matrix(NA, ncol = 2, nrow = length(all.genes))
	rownames(num.table) <- all.genes
	colnames(num.table) <- c("SA", "Gr")
	for(i in 1:length(all.genes)){
		num.table[i,1] <- mean(rank.table.SA[which(id.table.SA == all.genes[i])])
		num.table[i,2] <- mean(rank.table.Gr[which(id.table.Gr == all.genes[i])])
		}
	
	num.table[which(!is.finite(num.table))] <- 0
	plot(num.table[,1], num.table[,2], xlab = "Times Found by SA", ylab = "Times Found by Greedy")
	abline(0, 1)
	}