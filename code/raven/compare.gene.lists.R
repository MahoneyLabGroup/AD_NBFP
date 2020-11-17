#This function compares the ranked gene lists 


compare.gene.lists <- function(gene.stats1, gene.stats2, SA = TRUE){
	
	if(SA){
		table1 <- gene.stats1$id.table.SA
		table2 <- gene.stats2$id.table.SA		
		}else{
		table1 <- gene.stats1$id.table.Gr
		table2 <- gene.stats2$id.table.Gr		
		}
	
	n.loci <- nrow(table1)
	stat.table <- matrix(NA, nrow = n.loci+1, ncol = 6)
	colnames(stat.table) <- c("num.genes1", "num.genes2", "unique.genes", "Set1.not.in.Set2", "Set2.not.in.Set1", "Jaccard.Index")
	rownames(stat.table) <- c(rownames(table1), "overall")

	for(i in 1:(n.loci+1)){
		if(i > n.loci){
			gene.set1 <- table1[which(!is.na(table1))]
			gene.set2 <- table2[which(!is.na(table2))]			
			}else{
			gene.set1 <- table1[i,which(!is.na(table1[i,]))]
			gene.set2 <- table2[i,which(!is.na(table2[i,]))]
			}
		stat.table[i,1] <- length(gene.set1)
		stat.table[i,2] <- length(gene.set2)
		total.genes <- length(unique(c(gene.set1, gene.set2)))
		stat.table[i,3] <- total.genes
		shared.genes <- length(intersect(gene.set1, gene.set2))
		stat.table[i,4] <- length(setdiff(gene.set1, gene.set2))
		stat.table[i,5] <- length(setdiff(gene.set2, gene.set1))
		stat.table[i,6] <- shared.genes/total.genes
		}
		
	return(stat.table)
	
}