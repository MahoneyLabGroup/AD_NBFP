get.module.size <- function(results.dir, type = c("by.module", "by.set"), 
enrichment = FALSE, organism = "mmusculus", top.enrich = 5, verbose = FALSE){
	
	if(type == "by.set"){enrichment = FALSE}
	
	module.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.info[[1]]
	module.table <- module.info[[2]]
	
	num.genes <- rep(NA, nrow(module.table))
	enrich.terms <- rep(NA, nrow(module.table))
	for(i in 1:nrow(module.table)){
		if(verbose){cat(module.table[i,], "\n")}
		mod.genes <- read.csv(paste0(module.dir[i], "/Module.Gene.Info.csv"), stringsAsFactors = FALSE)
		if(enrichment){
			enrich <- gprofiler(mod.genes[,1], organism, src_filter = "GO")
			enrich.terms[i] <- paste(enrich[order(enrich[,"p.value"])[1:top.enrich],"term.name"], collapse = "; ")
			}
		num.genes[i] <- nrow(mod.genes)
		}
	
	if(enrichment){
		results <- cbind(module.table, num.genes, enrich.terms)
		colnames(results) <- c("Term", "Module", "N.Genes", "Enrichment.Terms")
		rownames(results) <- NULL
		}else{
		results <- cbind(module.table, num.genes)
		colnames(results) <- c("Term", "Module", "N.Genes")
		rownames(results) <- NULL		
		}
	
	if(type == "by.module"){
		return(results)
		}else{
		u_types <- unique(results[,1])
		Total <- rep(NA, length(u_types))
		num.types <- sapply(u_types, function(x) length(which(results[,1] == x)))
		new.table <- matrix("-", nrow = length(u_types), ncol = max(num.types))
		colnames(new.table) <- paste0("Module", 1:ncol(new.table))
		rownames(new.table) <- u_types
		for(i in 1:nrow(new.table)){
			type.locale <- which(results[,1] == u_types[i])
			mod.sizes <- results[type.locale,3]
			new.table[i,1:length(mod.sizes)] <- mod.sizes
			Total[i] <- sum(as.numeric(mod.sizes))
			}
		final.table <- cbind(Total, new.table)
		return(final.table)
		}
	

}