#This function calculates a false positive score for 
#all elements in a model based on their prediction scores
#predicted.decision.values should be a matrix with one column
#and rownames corresponding to the elements being predicted

svm.prediction.fp <- function(predicted.decision.values, svm.model){
	
	tp <- gsub("X", "", names(svm.model$x.scale[[1]]))
	all.genes <- names(svm.model$fitted)
	tn <- setdiff(all.genes, tp)
	
	decision.values <- svm.model$decision.values
	
	min.val <- min(decision.values[,1])
	max.val <- max(decision.values[,1])
	
	get.fp.tp <- function(test.val){
		pred.pos <- rownames(decision.values)[which(decision.values > test.val)]
		pred.neg <- rownames(decision.values)[which(decision.values <= test.val)]
		
		num.tp <- length(intersect(pred.pos, tp))
		num.tn <- length(intersect(pred.neg, tn))
		num.fp <- length(intersect(pred.pos, tn))
		num.fn <- length(intersect(pred.neg, tp))
		
		fp.rate <- num.fp/(num.fp + num.tn)
		tp.rate <- num.tp/(num.tp + num.fn)		
		
		result <- c("fp.rate" = fp.rate, "tp.rate" = tp.rate)
		return(result)
		}

	all.fp.tp <- lapply(predicted.decision.values[,1], get.fp.tp)
	fp.tp.mat <- Reduce("rbind", all.fp.tp)
	rownames(fp.tp.mat) <- rownames(predicted.decision.values)
		
	return(fp.tp.mat)
	
}