#This function is applied to a directory in which there is an 
#SVM.Prediction.Mat.RData file produced by training SVMs
#it calculates an ROC curve for all models in the path to indicate
#how well each SVM model separates the gene set from the rest
#of the genome
#The path specified as results.dir is the outer project results path. 
#It is assumed that in this directory, there are multiple gene sets
#and that each of these gene sets has one or more modules
#within it. This will be true if using generate.triage.models

svm.ROC <- function(results.dir = ".", plot.results = TRUE, n.cores = 4){

	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table	
			
	full.results <- vector(mode = "list", length = length(module.dir))
	names(full.results) <- module.dir


	for(i in 1:nrow(dir.table)){
		
		module.mat.file <- list.files(path = module.dir[i], pattern = "decomp", full.names = TRUE)
		if(length(module.mat.file) == 0){
			stop("generate.triage.models() must be run first to create the non.decomp.mat.RData 
			or the decomp.mat.RData file.")		
			}
			
		module.mat <- readRDS(module.mat.file)

		module.mat.base <- tail(unlist(strsplit(module.mat.file, "/")), 1)
		if(module.mat.base == "decomp.mat.RData"){
			#if an SVD was performed, read it in to determine the number of 
			#true positive genes
			decomp <- readRDS(paste0(module.dir[i], "/Decomposition.RData"))
			num.pos.genes <- length(decomp$d)
			}else{
			#otherwise, if no decomposition was done, the number of positives
			#is the number of rows in the module.mat
			num.pos.genes <- nrow(module.mat)	
			}
		
		#the positive genes are the first genes in the module matrix
		pos.genes <- colnames(module.mat)[1:num.pos.genes]
		neg.genes <- setdiff(colnames(module.mat), pos.genes)
	
		prediction.mat.file <- paste0(module.dir[i], "/SVM.Prediction.Mat.RData")
		if(!file.exists(prediction.mat.file)){
			stop("generate.triage.models() must be run first to create the 
			SVM.Prediction.Mat.RData file.")
			}
		prediction.mat <- readRDS(prediction.mat.file)
		
		prediction.fptp.file <- paste0(module.dir[i], "/Prediction.FP.TP.RData")
		if(!file.exists(prediction.fptp.file)){
			chunked.p <- chunkV(1:nrow(prediction.mat), n.cores)
			cl <- makeCluster(n.cores)
			registerDoParallel(cl)
			predict.fp.tp.list <- foreach(p = 1:length(chunked.p), .export = "svm.prediction.roc") %dopar% {
				#n.samples is the number of points along which to sample the ROC curve
				lapply(chunked.p[[p]], function(x) svm.prediction.roc(svm.prediction.scores = t(prediction.mat[x,,drop=FALSE]), true.positives = pos.genes, true.negatives = neg.genes, n.samples = 100, plot.results = FALSE))
				}				
			stopCluster(cl)
			predict.fp.tp <- unlist(predict.fp.tp.list, recursive = FALSE)
			saveRDS(predict.fp.tp, prediction.fptp.file)
			}else{
			predict.fp.tp <- readRDS(prediction.fptp.file)	
			}
			
		if(plot.results){
			smooth.auc(Reduce("rbind", predict.fp.tp), plot.results = TRUE)
			module.text <- paste(dir.table[i,], collapse = ": ")
			mtext(text = paste(module.text, "\n", length(pos.genes), "genes"))
			}
	
		full.results[[i]] <- predict.fp.tp
		}
		
		invisible(full.results)
	}