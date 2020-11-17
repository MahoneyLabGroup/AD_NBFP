#This function tells you if two segments overlap
segments.overlap <- function(start1, end1, start2, end2){	
	all.pos <- c("1" = start1, "1" = end1, "2" = start2, "2" = end2)
	sorted.pos <- sort(all.pos)
	does.overlap <- names(sorted.pos)[1] != names(sorted.pos)[2]
	return(does.overlap)
	}