#This function gets all the module directories in a project
#results directory (top.dir).
#if dir.table is TRUE, this function also returns a
#table of the last two directories in the hierarchy
#which included the gene list name and the module name

get.module.dir <- function(top.dir, dir.table = FALSE){
	
	gene.list.dir <- list.files(top.dir, full.names = TRUE)
	module.dir <- unlist(sapply(gene.list.dir, function(x) list.files(x, pattern = "Module", full.names = TRUE)))
	rdata.idx <- grep("RData", module.dir)
	if(length(rdata.idx) > 0){
		module.dir <- module.dir[-rdata.idx]
		}
	
	if(!dir.table){
		return(module.dir)
		}else{	
		split.dir <- strsplit(module.dir, "/")
		dir.table <- Reduce("rbind", lapply(split.dir, function(x) tail(x, 2)))
		if(is.null(dim(dir.table))){dir.table <- matrix(dir.table, nrow = 1)}
		final.results <- list("module.dir" = module.dir, "dir.table" = dir.table)
		return(final.results)
		}

	
	
}