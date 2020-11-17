#This function returns a lits of elements to use
#for a sliding window

sliding.window.el <- function(elV, window.size, gap.size){
	
	total.num = length(elV)
	el.list <- list()
	
	start.pos <- 1
	list.ind <- 1
	while(start.pos + window.size <= total.num){
		el.list[[list.ind]] <- elV[start.pos:(start.pos+window.size-1)]
		list.ind = list.ind + 1
		start.pos = start.pos + gap.size
		}
	if(start.pos < total.num){
		el.list[[list.ind]] <- elV[start.pos:total.num]
		}
	
	return(el.list)
	
	
	
	
}
