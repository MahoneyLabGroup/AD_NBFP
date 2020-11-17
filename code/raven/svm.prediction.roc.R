#This function takes in SVM results
#and plots an ROC curve
#true.positives <- mp.genes$entrez_id
#test.genes <- region.entrez
#you can select individual genes to calculate
#false positive scores for. For example,
#all genes in a given locus
#n.samples indicates how densely you want to sample the
#ROC curve. if n.samples = 1000, there will be 1000
#points on the curve sampled 
svm.prediction.roc <- function(svm.prediction.scores, true.positives, true.negatives, n.samples = 1000, plot.results = FALSE){
	
	test.val <- seq(min(svm.prediction.scores), max(svm.prediction.scores), length.out = n.samples)
	
	fp.rate <- rep(NA, length(test.val))
	tp.rate <- rep(NA, length(test.val))

	get.fp.tp <- function(test.val){
		pred.pos <- rownames(svm.prediction.scores)[which(svm.prediction.scores > test.val)]
		pred.neg <- rownames(svm.prediction.scores)[which(svm.prediction.scores <= test.val)]
		
		num.tp <- length(intersect(pred.pos, true.positives))
		num.tn <- length(intersect(pred.neg, true.negatives))
		num.fp <- length(intersect(pred.pos, true.negatives))
		num.fn <- length(intersect(pred.neg, true.positives))		
		
		fp.rate <- num.fp/(num.fp + num.tn)
		tp.rate <- num.tp/(num.tp + num.fn)		
		
		result <- list("fp.tp" = c("fp.rate" = fp.rate, "tp.rate" = tp.rate), c(num.tp, num.tn, num.fp, num.fn))
		return(result)
		}

	all.fp.tp <- lapply(test.val, get.fp.tp)
	fp.tp.mat <- Reduce("rbind", lapply(all.fp.tp, function(x) x[[1]]))
			
	if(plot.results){
		plot(fp.tp.mat[,1], fp.tp.mat[,2], pch = 16, ylab = "True Positive Rate", xlab = "False Positive Rate", xlim = c(0,1), ylim = c(0,1))
		abline(0, 1)
		}
	
	return(fp.tp.mat)
	
}