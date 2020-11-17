#Retrieves the SVM score matrices from all modules in a results
#directory. If reduce == TRUE, the matrices will be reduced to only
#the common genes between them. Otherwise, the full matrices are returned

get.prediction.mats <- function(results.dir, reduce = FALSE){

	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table	

	all.predictions <- vector(mode = "list", length = length(module.dir))
	names(all.predictions) <- apply(dir.table, 1, function(x) paste(x, collapse = "_"))
	for(i in 1:length(module.dir)){
		pred.mat.file <- paste0(module.dir[i], "/SVM.Prediction.Mat.RData")	
		all.predictions[[i]] <- readRDS(pred.mat.file)
		}

	get.common.genes <- function(gene.list, pred.mat){
		gene.locale <- match(gene.list, colnames(pred.mat))
		red.mat <- pred.mat[,gene.locale]
		return(red.mat)
		}

	if(reduce){
		common.genes <- Reduce("intersect", lapply(all.predictions, colnames))
		red.mats <- lapply(all.predictions, function(x) get.common.genes(common.genes, x))
		}else{
		red.mats <- all.predictions	
		}


	return(red.mats)

	}