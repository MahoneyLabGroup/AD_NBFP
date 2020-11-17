#This function takes a tissue matrix and performs SVM
#vote.min is the minimum number of positive votes required
#for a gene to be considered a postive
#one of prediction.mat or trial.num needs to be specified
#prediction.mat is the output of tissue.svm, but this
#function can also read the SVM model ouptut from that
#function using trial.num

predict.genes.in.class <- function(path = ".", tissue.mat, prediction.mat = NULL, trial.num = NULL, vote.min = 0.5){
	
	prediction.file = "SVM.Prediction.Mat.RData"
	
	if(is.null(prediction.mat) && is.null(trial.num)){stop("One of prediction.mat or trial.num must be specified.\n")}
	if(!is.null(prediction.mat) && !is.null(trial.num)){stop("Only one of prediction.mat or trial.num may be specified.\n")}


	tissue.mat <- as.matrix(tissue.mat)
	
	tp.mat <- tissue.mat[,1:nrow(tissue.mat)]
	tn.mat <- tissue.mat[,((nrow(tissue.mat)+1):ncol(tissue.mat))]
	
	full.data <- cbind(tp.mat, tn.mat)
	labels <- as.factor(c(rep("pos", ncol(tp.mat)), rep("neg", ncol(tn.mat))))

	num.tp <- ncol(tp.mat)
	num.tn <- ncol(tn.mat)
	
	sample.labels <- as.factor(c(rep("pos", num.tp), rep("neg", num.tp)))

	if(!is.null(trial.num)){
		all.predict <- matrix(NA, nrow = length(trial.num), ncol = ncol(full.data))
		colnames(all.predict) <- colnames(full.data)
		
		for(n in 1:length(trial.num)){
			svm.results <- readRDS(paste0(path, "/SVM_Results", trial.num[n], ".RData"))		
			
			prediction <- predict(svm.results$opt.model, newdata = t(full.data), decision.values = TRUE)
			prediction.order <- match(rownames(attributes(prediction)$"decision.values"), colnames(all.predict))
			
			all.predict[n,] <- as.vector(attributes(prediction)$"decision.values"[prediction.order])
			}
		saveRDS(all.predict, paste0(path, "/", prediction.file))
		}else{
		all.predict <- prediction.mat
		}

	
}