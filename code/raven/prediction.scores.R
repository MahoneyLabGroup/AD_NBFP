#This function retrieves the mean prediction
#scores for all genes in all modules in a results 
#directory

prediction.scores <- function(results.dir = "."){

	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table	
		
	full.results <- vector(mode = "list", length = length(module.dir))
	names(full.results) <- apply(dir.table, 1, function(x) paste(x, collapse = "_"))


	for(i in 1:nrow(dir.table)){
			
		prediction.mat.file <- paste0(module.dir[i], "/SVM.Prediction.Mat.RData")
		if(!file.exists(prediction.mat.file)){
			stop("generate.triage.models() must be run first to create the SVM.Prediction.Mat.RData file.")
			}
		prediction.mat <- readRDS(prediction.mat.file)
	
		full.results[[i]] <- colMeans(prediction.mat)
		}
		
		invisible(full.results)
	}