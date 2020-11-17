#This function takes a tissue matrix and performs SVM


tissue.svm <- function(path = ".", tissue.mat, num.tp = nrow(tissue.mat), trial.num = 1:100, C.list = 4^seq(-5, -3, 1), verbose = FALSE){
	
	tissue.mat <- as.matrix(tissue.mat)
	
	tp.mat <- tissue.mat[,1:num.tp]
	tn.mat <- tissue.mat[,(num.tp+1):ncol(tissue.mat)]
	
	full.data <- cbind(tp.mat, tn.mat)
	labels <- as.factor(c(rep("pos", ncol(tp.mat)), rep("neg", ncol(tn.mat))))

	num.tp <- ncol(tp.mat)
	num.tn <- ncol(tn.mat)
	
	all.predict <- matrix(NA, nrow = length(trial.num), ncol = ncol(full.data))
	colnames(all.predict) <- colnames(full.data)
	rownames(all.predict) <- trial.num

	sample.labels <- as.factor(c(rep("pos", num.tp), rep("neg", num.tp)))

	for(n in 1:length(trial.num)){
		cat("Trial", trial.num[n], "\n")
		sample.mat <- t(cbind(tp.mat, tn.mat[,sample(1:num.tn, num.tp)]))
		
		results <- cv.linear.svm(data.mat = sample.mat, data.labels = sample.labels, C.list = C.list, verbose = verbose)
		
		if(is.null(results$opt.model)){
			return(NULL)
			}
		
		prediction <- predict(results$opt.model, newdata = t(full.data), decision.values = TRUE)
		prediction.order <- match(rownames(attributes(prediction)$"decision.values"), colnames(all.predict))
		
		all.predict[n,] <- as.vector(attributes(prediction)$"decision.values"[prediction.order])

		saveRDS(results, paste0(path, "/SVM_Results", trial.num[n], ".RData"))
		}
	
	return(all.predict)
	
}