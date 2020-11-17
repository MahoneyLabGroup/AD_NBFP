#I am starting to write a function to get genes associated with
#MP terms using ontoCAT. The code here so far is from ontology.testing.R
#THIS NEEDS UPDATING FOR HUMANS AND YEAST TO MAKE A TABLE FOR THE RESULTS INSTEAD OF A LIST


get.mp.genes <- function(mp.term, organism = c("mouse", "human", "yeast"), verbose = FALSE){

	require(InterMineR)
	require(biomaRt)
	organism <- organism[1]

	all.var <- ls(globalenv())
	
	if(organism[1] == "mouse"){	
		mine <- initInterMine("http://www.mousemine.org/mousemine")
		# getTemplates(mine)
		gene.query <- getTemplateQuery(mine, "MPhenotype_MGenotype")
		gene.query$select <- c(gene.query$select, "OntologyAnnotation.subject.primaryIdentifier")
		value.idx <- 5

		if(length(which(all.var == "lib")) == 0){
			# lib <<- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
			lib <<- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
			}

		}
		
	if(organism[1] == "human"){
		mine <- initInterMine("http://www.humanmine.org/humanmine")
		gene.query <- getTemplateQuery(mine, "PhenotypeGene")
		value.idx <- 2

		if(length(which(all.var == "lib")) == 0){
		# lib <<- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
		lib <<- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
			}		
		}
	
	if(organism[1] == "yeast"){
		mine <- initInterMine("https://yeastmine.yeastgenome.org/yeastmine")
		gene.query <- getTemplateQuery(mine, "Phenotype_Genes")
		value.idx <- 1
		}

	results <- NULL
	for(i in 1:length(mp.term)){
		if(verbose){cat("Looking up", mp.term[i], "...\n")}
		if(mp.term[i] == "*"){
			gene.query$where[[value.idx]]$value <- "*"
			}else{
			gene.query$where[[value.idx]]$value <- paste0("*", mp.term[i], "*")
			}
		results <- rbind(results, runQuery(mine, gene.query))
		}
	if(is.null(results)){return("none")}
	
	if(organism == "mouse"){
		gene.names <- unique(results[,"OntologyAnnotation.subject.symbol"])
		allele.locale <- grep("\\(", gene.names)
		gene.names = gene.names[setdiff(1:length(gene.names), allele.locale)]
		
		if(length(gene.names) == 0){return("none")}
		gene.table <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", gene.names, lib)
		no.nas <- gene.table[which(!is.na(gene.table[,2])),]
		gene.mp <- lapply(no.nas[,1], function(x) results[which(results[,"OntologyAnnotation.subject.symbol"] == x),"OntologyAnnotation.ontologyTerm.name"])
		MP.term <- unlist(lapply(gene.mp, function(x) paste(unique(x), collapse = ";")))
		final.table <- cbind(no.nas, MP.term)
		final.result <- final.table[order(final.table[,1]),]
		}
		
	if(organism == "human"){
		gene.table <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", results[,"Gene.homologues.homologue.symbol"], lib)
		not.na <- which(!is.na(gene.table[,2]))
		final.result <- list("terms" = unique(results[,"Gene.alleles.genotypes.phenotypeTerms.name"]), "gene.names" = unique(gene.table[not.na,1]), "gene.id" = unique(gene.table[not.na,2]))
		}
		
		
	if(organism == "yeast"){
		gene.name <- results[,"Phenotype.genes.symbol"]
		gene.id <- results[,"Phenotype.genes.secondaryIdentifier"]	
		terms <- results[ ,"Phenotype.observable"]
		final.result <- list("terms" = unique(terms), "gene.names" = unique(gene.name), "gene.id" = unique(gene.id))
		}
	
	if(verbose){
	if(class(final.result) == "list"){
		cat("I found\n")
		cat(length(final.result[[1]]), "unique terms and\n")
		cat(length(final.result[[3]]), "unique gene IDs\n")	
		}else{
		cat("I found\n")
		cat(length(unique(unlist(gene.mp))), "unique terms and\n")
		cat(nrow(final.result), "unique gene IDs\n")				
		}
	}
	
	return(final.result)
	
	}