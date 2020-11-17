#This function tunes the cost parameter by tracking cross validation accuracy
#over a range of costs
#It is now required to put in the number of true positive genes
#since sometimes we will be using eigenvectors instead of individual genes

tune.cost.parameter <- function(tissue.adj, num.tp = nrow(tissue.adj), n.trials = 10, C.list = 10^seq(-5, 1, 1), verbose = TRUE){
	
	tissue.adj <- as.matrix(tissue.adj)

	tp.mat <- tissue.adj[,1:num.tp]
	tn.mat <- tissue.adj[,(num.tp+1):ncol(tissue.adj)]
	
	full.data <- cbind(tp.mat, tn.mat)
	labels <- as.factor(c(rep("pos", ncol(tp.mat)), rep("neg", ncol(tn.mat))))

	num.tp <- ncol(tp.mat)
	num.tn <- ncol(tn.mat)
	
	all.accuracy <- matrix(0, nrow = n.trials, ncol = length(C.list))
	#generate balanced sample labels
	sample.labels <- as.factor(c(rep("pos", num.tp), rep("neg",num.tp)))
	
	acc.means <- colMeans(all.accuracy)
	old.max <- max(acc.means)
	new.max <- 10
	run.num <- 1
	new.C.list <- C.list
	
	while(new.max > old.max){
		changed.C <- 0
		old.max <- new.max
			
		for(n in 1:n.trials){
			if(verbose){
				cat("\nTrial", n, "\n")
				}
			sample.mat <- t(cbind(tp.mat, tn.mat[,sample(1:num.tn, num.tp)]))
			
			results <- cv.linear.svm(data.mat = sample.mat, data.labels = sample.labels, C.list = new.C.list, verbose = verbose)
			all.accuracy[n,] <- results[[2]]
			}
				
		acc.means <- colMeans(all.accuracy)
		new.max <- max(acc.means)
		max.acc <- which.max(acc.means)
		
		#create a new cost list based on where the maximum accuracy is
		#going over the initial cost range is causing problems because
		#higher costs take much more time to compute. I'm now restricting
		#the cost list to the limits of the original list
		
		#if the maximum accuracy is at the bottom of the list
		if(max.acc == 1){
			min.c <- new.C.list[max.acc]/10
			max.c <- new.C.list[(max.acc+1)]
			changed.C <- 1
			}
		#if the maximum accuracy is at the top of the list
		if(max.acc == length(new.C.list)){
			min.c <- new.C.list[(max.acc-1)]
			max.c <- new.C.list[max.acc]
			changed.C <- 1
			}
		if(changed.C == 0){
			min.c <- new.C.list[(max.acc-1)]
			max.c <- new.C.list[(max.acc+1)]
			}				
		#generate a new sequence between the new min and max
		new.C.list <- segment.region(min.c, max.c, length(new.C.list), "ends")
		}
		
	return(new.C.list)	
		
	}