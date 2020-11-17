#This function runs Hotelling's T test for a phenotype 
#group on a per-marker basis using the R package Hotelling.

marker.hotelling <- function(geno.mat, pheno.mat, verbose = FALSE){
	require(Hotelling)
		
	all.tests <- vector(mode = "list", length = ncol(geno.mat))
	x <- data.frame(pheno.mat)
	for(i in 1:ncol(geno.mat)){
		if(verbose){report.progress(i, ncol(geno.mat))}
		y <- data.frame(geno.mat[,i,drop=FALSE])
		colnames(y) <- "y"
		df <- cbind(y,x)
		all.tests[[i]] <- hotelling.test(.~y, data = df, shrinkage = TRUE) 	
	}
	return(all.tests)
}
