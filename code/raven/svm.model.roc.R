svm.model.roc <- function(svm.model, perm = 100, plot.results = FALSE){
	
	if(!plot.results){
		perm = 0
		}
	
	
	n.genes <- nrow(svm.model$decision.values)
	tp <- rownames(svm.model$decision.values)[1:(n.genes/2)]
	# tp <- gsub("X", "", names(svm.model$x.scale[[1]]))
	
	all.genes <- names(svm.model$fitted)
	tn <- setdiff(all.genes, tp)
	
	decision.values <- svm.model$decision.values
	
	min.val <- min(decision.values[,1])
	max.val <- max(decision.values[,1])
	
	val.seq <- seq(min.val, max.val, 0.1)
	
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

	all.fp.tp <- lapply(val.seq, get.fp.tp)
	fp.tp.mat <- Reduce("rbind", all.fp.tp)
	
	if(plot.results){
		plot(fp.tp.mat[,1], fp.tp.mat[,2], pch = 16, ylab = "True Positive Rate", xlab = "False Positive Rate")
		abline(0, 1)
		}
	
	if(perm > 0){
		for(p in 1:perm){
			gene.perm <- sample(all.genes)
			tp <- gene.perm[1:(length(gene.perm)/2)]
			tn <- setdiff(gene.perm, tp)
			all.fp.tp <- lapply(val.seq, get.fp.tp)
			fp.tp.mat <- Reduce("rbind", all.fp.tp)
			if(plot.results){
				plot(fp.tp.mat[,1], fp.tp.mat[,2], pch = 16, ylab = "True Positive Rate", xlab = "False Positive Rate")
				abline(0, 1)
				}
			
			}

		}
	
	invisible(fp.tp.mat)
	
}