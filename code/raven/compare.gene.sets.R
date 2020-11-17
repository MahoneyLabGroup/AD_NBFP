#This function is used for comparing gene sets with with overlapping
#membership
#gene.list is a list containing the vectors to be compared. It doesn't
#matter what is in these vectors as long as the naming conventions
#are consistent across all elements of the list
#This function makes a barplot showing the number of genes in each
#list, a barplot showing the overlap of genes between each list pair
#and a heatmap of the jaccard indices of all list pairs


compare.gene.sets <- function(gene.list, plot.type = c("gene.num", "shared.genes", "jaccard")){
	
	#layout(matrix(c(1,2,0,2), nrow = 2, byrow = TRUE), widths = c(2,1.5))
	
	cols <- c("#d8b365", "#f5f5f5", "#5ab4ac")
	
	num.genes <- sapply(gene.list, length)
	num.order <- order(num.genes, decreasing = TRUE)
	
	if(length(grep("gene.num", plot.type)) == 1){
		par(mar = c(10, 4, 4, 4))
		barplot(num.genes[num.order], las = 2)
		}
	
	list.pairs <- pair.matrix(1:length(gene.list))
	shared.mat <- matrix(NA, ncol = 3, nrow = nrow(list.pairs))
	colnames(shared.mat) <- c("only1", "shared", "only2")
	total.genes <- rep(NA, nrow(list.pairs))
	jaccard.pairs <- matrix(NA, nrow = length(gene.list), ncol = length(gene.list))
	rownames(jaccard.pairs) <- colnames(jaccard.pairs) <- names(gene.list)
	for(i in 1:nrow(list.pairs)){
		set1 <- list.pairs[i,1]
		set2 <- list.pairs[i,2]
		shared.mat[i,1] <- length(setdiff(gene.list[[set1]], gene.list[[set2]]))
		shared.mat[i,2] <- length(intersect(gene.list[[set1]], gene.list[[set2]]))
		shared.mat[i,3] <- length(setdiff(gene.list[[set2]], gene.list[[set1]]))
		jaccard.pairs[set1, set2] <- jaccard.pairs[set2, set1] <- jaccard.ind(gene.list[[set1]], gene.list[[set2]])
		total.genes[i] <- length(unique(c(gene.list[[set1]], gene.list[[set2]])))
		}

	prop.mat <- apply(shared.mat, 2, function(x) x/total.genes)	

	comp.names <- apply(list.pairs, 1, function(x) paste(c(names(gene.list)[x[1]], names(gene.list)[x[2]]), collapse = "_"))
	shared.order <- order(prop.mat[,"shared"], decreasing = FALSE)
	par(mar = c(4, 20, 4, 4))

	if(length(grep("shared.genes", plot.type)) == 1){
	barplot(t(prop.mat[shared.order,]), col = cols, names = comp.names[shared.order], las = 2, horiz = TRUE)	
	}

	if(length(grep("jaccard", plot.type)) == 1){
	 pheatmap(jaccard.pairs)
	}

}