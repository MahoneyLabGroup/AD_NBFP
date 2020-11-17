#This function splits a set of string by a series of patterns

multi.strsplit <- function(strV, patterns){
	
	new.string <- strV
	for(i in 1:length(patterns)){
		split.strings <- strsplit(new.string, patterns[i])
		new.string <- unlist(lapply(split.strings, function(x) x[which(x != "")]))
		}
	
	return(new.string)
	
}