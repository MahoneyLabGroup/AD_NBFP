#This function compares two ranked lists
#It ranks elements in the first list from 1:N
#for the second list, it assigns the ranks from
#the first list. It then uses a Spearman correlation
#to compare the lists.
#this might take a while for very long lists. To only 
#compare the top of ranked.list2 to ranked, set top.N
#to the maximum rank of ranked.list2 that should be ranked
#relative to ranked.list1.

#test <-readRDS('~/Documents/Projects/FNTM_and_GIANT/results/SVM_mouse/DLPFCbrown/Module1/SVM.Prediction.Mat.RData') 
#ranked.list1 <- colnames(test)[order(colMeans(test[1:50,]))]
#ranked.list2 <- colnames(test)[order(colMeans(test[51:100,]))]

compare.ranked.lists <- function(ranked.list1, ranked.list2, top.N = NULL, plot.results = FALSE){

	if(is.null(top.N)){top.N = min(c(length(ranked.list2), length(ranked.list1)))}

	ranks1 <- 1:length(ranked.list1)

	get.rank <- function(value, the.list){
		val.locale <- which(the.list == value)
		if(length(val.locale) == 0){
			return(NA)
			}else{
			return(val.locale)
			}
	}

	ranks2 <- sapply(ranked.list2[1:top.N], function(x) get.rank(x, ranked.list1))	
	
	if(plot.results){
		plot(ranks1[1:top.N], ranks2[1:top.N], xlab = "Ranked List 1", ylab = "Ranked List 2")
		abline(0,1)
		}

	s.cor <- cor(ranks1[1:top.N], ranks2[1:top.N], method = "spearman", use = "pairwise.complete")	
	results <- list("ranks1" = ranks1, "ranks2" = ranks2, "spearman.cor" = s.cor)
	return(results)
		
}