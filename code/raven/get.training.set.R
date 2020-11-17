#This function retrieves all the true positive genes 
#in the submodules in a base directory.


get.training.set <- function(results.dir){

	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table	
    mod.names <- apply(dir.table, 1, function(x) paste(x[1], x[2], sep = "_"))

    all.tp <- lapply(module.dir, function(x) rownames(readRDS(file.path(x, "non.decomp.mat.RData"))))
    names(all.tp) <- mod.names
    return(all.tp)
}