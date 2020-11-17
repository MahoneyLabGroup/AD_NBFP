#This function takes in a gene and gets all relevant MP terms

get.mp.terms <- function(gene.name, organism = c("mouse", "human", "yeast")){

	
	require(ontoCAT)
	require(biomaRt)

	base.dir <- "~/Documents/Data/Little_Cross/little_cross_data_cape/Body_Comp_IGF_New_P/IMP"
	all.var <- ls(globalenv())


	organism <- organism[1]

	if(organism[1] == "mouse"){		
		obo.file = "file:/Users/atyler/Documents/Ontologies/MPheno_OBO.ontology.txt"
		mp.gene.file = "~/Documents/Ontologies/MGI_GenePheno.rpt.txt"
		term.col <- 5; gene.col <- 1
		if(length(which(all.var == "mp.mart")) == 0){
			mp.mart <<- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
			}
		}
		
	if(organism[1] == "human"){
		obo.file <- "file:/Users/atyler/Documents/Ontologies/hp.obo"
		mp.gene.file <- "~/Documents/Ontologies/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt"
		term.col <- 4; gene.col <- 2
		if(length(which(all.var == "mp.mart")) == 0){
			mp.mart <<- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
			}
		}
	
	if(organism[1] == "yeast"){
		obo.file <- "file:/Users/atyler/Documents/Ontologies/goslim_yeast.obo"
		mp.gene.file = "~/Documents/Ontologies/go_slim_mapping.txt"
		term.col <- 6; gene.col <- 2
		if(length(which(all.var == "mp.mart")) == 0){
			mp.mart <<- useEnsembl(biomart="ensembl", dataset="scerevisiae_gene_ensembl") 		
			}
		}
	
	
	#get phenotype ontology
	if(length(which(all.var == "mp.onto")) == 0){
		mp.onto <<- getOntology(obo.file)
		}
	
	#read in gene/pheno associations
	if(length(which(all.var == "geno.pheno")) == 0){
		gene.pheno <<- as.matrix(read.table("~/Documents/Ontologies/MGI_PhenotypicAllele.rpt.txt", sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
		}

	gene.locale <- grep(gene.name, gene.pheno[,8], ignore.case = TRUE)
	if(length(gene.locale) == 0){return("no terms found")}

	all.terms <- unlist(lapply(gene.pheno[gene.locale,11], function(x) strsplit(x, ",")))
	
	all.descriptions <- unlist(lapply(all.terms, function(x) getTermNameById(mp.onto, x)))
		
	return(all.descriptions)

}