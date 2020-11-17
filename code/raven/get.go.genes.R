#I am starting to write a function to get genes associated with
#MP terms using ontoCAT. The code here so far is from ontology.testing.R
#THIS NEEDS UPDATING FOR HUMANS AND YEAST TO MAKE A TABLE FOR THE RESULTS INSTEAD OF A LIST


get.go.genes <- function(go.term, organism = c("mouse", "human", "yeast")){


	org.options <- c("mouse", "human", "yeast")
	require(InterMineR)
	require(biomaRt)
	organism <- organism[1]

	org.locale <- which(org.options == organism)
	if(length(org.locale) == 0){
		cat("The possible organisms are mouse, yeast, or human. Please check spelling.\n")
		}

	all.var <- ls(globalenv())
	
	
	#========================================================
	# for mice
	#========================================================
	
	if(organism[1] == "mouse"){	
		mine <- initInterMine("http://www.mousemine.org/mousemine")
		# getTemplates(mine)
		
		#find all GO terms with string we are looking for
		term.search <- getTemplateQuery(mine, "Lookup_GO")
		term.search$where[[2]]$value <- go.term
		term.results <- runQuery(mine, term.search)
		all.go.terms <- term.results[,1]
		all.go.names <- term.results[,2]
			
		term.query <- getTemplateQuery(mine, "Term_Descendants")
		more.go.terms <- NULL
		more.go.names <- NULL
		for(i in 1:length(all.go.terms)){
			term.query$where[[1]]$value <- all.go.terms[i]
			more.terms <- runQuery(mine, term.query)
			more.go.terms <- c(more.go.terms, more.terms[,3])
			more.go.names <- c(more.go.names, more.terms[,4])
			}
		all.go.terms <- unique(c(all.go.terms, more.go.terms))
		all.go.names <- unique(c(all.go.names, more.go.names))
		
		if(length(all.go.terms) == 0){
			cat("I did not find any GO terms matching the string.\n")
			return(NULL)
			}else{
			cat("I am looking for genes associated with", length(all.go.terms), "GO terms.\n")				
			}

		
		#now get all genes associated with GO terms
		gene.query <- getTemplateQuery(mine, "GO_MFeatures_NoChildren")
		results <- NULL
		for(i in 1:length(all.go.terms)){
			gene.query$where[[4]]$value <- all.go.terms[i]
			results <- rbind(results, runQuery(mine, gene.query))
			}

		all.terms.found <- unique(results[,2])
		all.names.found <- unique(results[,3])
		all.genes.found <- unique(results[,8])
		
		if(length(all.genes.found) == 0){
			cat("I did not find any genes associated with these terms.\n")
			return(NULL)
			}
	
		
		if(length(which(all.var == "lib")) == 0){
			# lib <<- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
			lib <<- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
			}

		gene.table <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", all.genes.found, lib)
		no.nas <- gene.table[which(!is.na(gene.table[,2])),]
		
		#figure out how to make a table that groups the GO terms by gene. 
		#this is still the old code, and it doesn't work
		
		gene.go <- lapply(no.nas[,1], function(x) results[which(results[,"SequenceFeature.symbol"] == x),"SequenceFeature.ontologyAnnotations.ontologyTerm.name"])
		gene.term.text <- unlist(lapply(gene.go, function(x) paste(unique(x), collapse = ";")))
		
		gene.order <- match(all.genes.found, gene.table[,1])
		final.gene.table <- cbind(gene.table[gene.order,], gene.term.text)

		}



	#========================================================
	# for humans
	#========================================================
	if(organism[1] == "human"){
		mine <- initInterMine("http://www.humanmine.org/humanmine")
		# getTemplates(mine)
		
		#find all GO terms with string we are looking for
		go.term.query <- getTemplateQuery(mine, "GOname_GOidentifier")
		go.term.query$where[[1]]$value <- paste0("*", go.term, "*")
		term.results <- runQuery(mine, go.term.query)
		all.go.terms <- term.results[,2]
		all.go.names <- term.results[,1]
		
		if(length(all.go.terms) == 0){
			cat("I did not find any GO terms matching the string.\n")
			return(NULL)
			}else{
			cat("I am looking for genes associated with", length(all.go.terms), "GO terms.\n")				
			}
		
		#find genes associated with GO terms
		gene.query <- getTemplateQuery(mine, "GOterm_Gene")
		results <- NULL
		for(i in 1:length(all.go.terms)){
			gene.query$where[[1]]$value <- all.go.names[i]
			results <- rbind(results, runQuery(mine, gene.query))
			}

		all.terms.found <- unique(results[,4])
		all.names.found <- unique(results[,5])
		all.genes.found <- unique(results[,2])

		if(length(all.genes.found) == 0){
			cat("I did not find any genes associated with these terms.\n")
			return(NULL)
			}
	
		if(length(which(all.var == "lib")) == 0){
		# lib <<- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
		lib <<- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
			}		
		
		gene.table <- getBM(c("external_gene_name", "entrezgene"), "external_gene_name", all.genes.found, lib)
		not.na <- which(!is.na(gene.table[,2]))
		gene.terms <- lapply(all.genes.found, function(x) unique(results[which(results[,"Gene.symbol"] == x),"Gene.goAnnotation.ontologyTerm.name"]))
		gene.term.text <- unlist(lapply(gene.terms, function(x) paste(x, collapse = ";")))
		gene.order <- match(all.genes.found, gene.table[,1])
		final.gene.table <- cbind(gene.table[gene.order,], gene.term.text)
		
		}

	#========================================================
	# for yeast
	#========================================================
	
	if(organism[1] == "yeast"){
		mine <- initInterMine("https://yeastmine.yeastgenome.org/yeastmine")
		# getTemplates(mine)
		
		#find all GO terms with string we are looking for
		go.term.query <- getTemplateQuery(mine, "GOname_GOidentifier")
		go.term.query$where[[3]]$value <- paste0("*", go.term, "*")
		term.results <- runQuery(mine, go.term.query)
		all.go.terms <- term.results[,2]
		all.go.names <- term.results[,1]

		if(length(all.go.terms) == 0){
			cat("I did not find any GO terms matching the string.\n")
			return(NULL)
			}else{
			cat("I am looking for genes associated with", length(all.go.terms), "GO terms.\n")				
			}

		#find genes associated with GO terms
		gene.query <- getTemplateQuery(mine, "GOTerm_Genes")
		results <- NULL
		for(i in 1:length(all.go.names)){
			gene.query$where[[4]]$value <- all.go.names[i]
			results <- rbind(results, runQuery(mine, gene.query))
			}

		all.terms.found <- unique(results[,5])
		all.names.found <- unique(results[,6])
		all.genes.found <- unique(results[,2])

		if(length(all.genes.found) == 0){
			cat("I did not find any genes associated with these terms.\n")
			return(NULL)
			}
		
		gene.terms <- lapply(all.genes.found, function(x) unique(results[which(results[,"Gene.secondaryIdentifier"] == x),"Gene.goAnnotation.ontologyTerm.name"]))
		gene.term.text <- unlist(lapply(gene.terms, function(x) paste(x, collapse = ";")))
		gene.names <- unlist(lapply(all.genes.found, function(x) unique(results[which(results[,"Gene.secondaryIdentifier"] == x),"Gene.symbol"])))

		final.gene.table <- cbind(gene.names, all.genes.found, gene.term.text)	
		}
	
	colnames(final.gene.table) <- c("gene.name", "secondary.identifier", "associated.GO.terms")
	
	cat("I found\n")
	cat(length(all.terms.found), "unique terms and\n")
	cat(length(all.genes.found), "unique gene IDs\n")				
	
	return(final.gene.table)
	
	}