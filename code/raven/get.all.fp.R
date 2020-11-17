get.all.fp <- function(results.dir){
	
	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table
		
	all.fp <- vector(mode = "list", length = length(module.dir))
	
	for(i in 1:length(module.dir)){
		results.file <- paste0(module.dir[i], "/Candidate.Gene.Results.csv")
		results <- read.csv(results.file, stringsAsFactors = FALSE)
	
		fp <- results[,c("external_gene_name","Mean.FP.Rate")]
		colnames(fp)[2] <- paste(dir.table[i,], collapse = "_")
		fp[,1] <- rename.dups(fp[,1], "character")
		
		all.fp[[i]] <- fp
	
		}

	all.gene.names <- unique(unlist(lapply(all.fp, function(x) x[,1])))
	
	fp.mat <- matrix(NA, nrow = length(all.gene.names), ncol = length(module.dir))
	rownames(fp.mat) <- all.gene.names
	colnames(fp.mat) <- 1:ncol(fp.mat)
	for(i in 1:length(module.dir)){
		gene.locale <- match(all.fp[[i]][,1], all.gene.names)
		fp.mat[gene.locale,i] <- all.fp[[i]][,2]
		colnames(fp.mat)[i] <- colnames(all.fp[[i]])[2]
		}

	return(fp.mat)
	
	}