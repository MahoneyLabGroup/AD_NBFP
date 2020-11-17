#This function builds the constraint matrix for 
#the quadratic programming package quadprog

constraint.matrices <- function(locus.genes){
	
	#The equality matrix tells us which genes are in which loci
	all.genes <- unlist(locus.genes)
	equality.matrix <- matrix(0, ncol = length(all.genes), nrow = length(locus.genes))
	rownames(equality.matrix) <- names(locus.genes)
	colnames(equality.matrix) <- all.genes
	
	for(i in 1:length(locus.genes)){
		equality.matrix[i,as.character(locus.genes[[i]])] <- 1
		}
	
	equality.constraint <- rep(1, length(locus.genes))
	
	
	#The first inequality matrix says we want all values of x to be greater or equal to 0
	inequality.matrix1 <- diag(1, length(all.genes), length(all.genes))
	inequality.v1 <- rep(0, length(all.genes))
	
	inequality.matrix2 <- diag(-1, length(all.genes), length(all.genes))
	inequality.v2 <- rep(-1, length(all.genes))
	
	Amat <- rbind(equality.matrix, inequality.matrix1, inequality.matrix2)
	bvec <- matrix(c(equality.constraint, inequality.v1, inequality.v2), ncol = 1)
	
	final.result <- list("Amat" = t(Amat), "bvec" = bvec)
	
}