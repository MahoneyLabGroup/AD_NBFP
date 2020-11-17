#This function plots a precision recall curve


pc.curve <- function(true.class, num.votes){
	
	max.votes <- max(num.votes)
	tp <- which(true.class == "pos")
	tn <- which(true.class == "neg")

	precision <- rep(NA, max.votes)
	recall <- rep(NA, max.votes)

	for(i in 1:max.votes){
		above.cutoff <- which(num.votes >= i)
		below.cutoff <- which(num.votes < i)
		
		tp.i <- intersect(tp, above.cutoff)
		tn.i <- intersect(tn, below.cutoff)
		
		fp.i <- intersect(tn, above.cutoff)
		fn.i <- intersect(tp, below.cutoff)
		
		precision[i] <- length(tp.i)/(length(tp.i)+length(fp.i))
		recall[i] <- length(tp.i)/(length(tp.i)+length(fn.i))
		
		}

	recall[which(!is.finite(recall))] <- 0
	precision[which(!is.finite(precision))] <- 0
	
	plot(recall, precision, type = "l", xlim = c(0, 1), ylim = c(0,1))

	# plot(precision, type = "l", ylim = c(0,1))
	# plot(recall, type = "l", ylim = c(0,1))
	# cbind(recall, precision)

}