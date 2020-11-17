#This function cycles through phenotype search
#terms and searches the MPO database to build
#an MP term - gene object

build.mp.obj <- function(phenotype.terms, organism = c("mouse", "human", "yeast")){
	
	organism <- organism[1]

	all.mp.obj <- vector(mode = "list", length = length(phenotype.terms))
	names(all.mp.obj) <- phenotype.terms
	
	for(ph in 1:length(phenotype.terms)){
		cat("\n", ph, phenotype.terms[ph], "\n")
		mp.genes <- get.mp.genes(phenotype.terms[ph], organism, "entrezgene")
		all.mp.obj[[ph]] <- mp.genes
		}

	return(all.mp.obj)	
	
}