#This function reads in specific lines from a file
#This is useful in situations in which a file is 
#too large to store in memory.


# filename <- "rankMatrix.txt"; line.num = c(4,6,9,10)
read.lines <- function(filename, line.num, delim = "\t", na.strings = "NA", blank.lines.skip = TRUE){
	
	the.lines <- apply(matrix(line.num, ncol = 1), 1, function(x) read.table(filename, sep = delim, na.strings = na.strings, skip = (x-1), nrows = 1, stringsAsFactors = FALSE, blank.lines.skip = blank.lines.skip))
	mat <- as.matrix(t(sapply(the.lines, function(x) as.vector(x))))
	
	# start.time <- Sys.time()
	# the.lines <- t(apply(matrix(line.num, ncol = 1), 1, function(x) scan(filename, sep = delim, na.strings = na.strings, skip = (x-1), nlines = 1, what = "character")))
	# end.time <- Sys.time()
	# total.time <- end.time - start.time
	
	return(mat)
	
	
}