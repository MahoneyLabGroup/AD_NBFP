#This function returns a number with padding
#0's to make the number a fixed number of digits

fixed.digit <- function(x, digits = 2){
	
	num.digs <- length(strsplit(as.character(x), "")[[1]])
	if(num.digs > digits){
		stop("The value entered alread exceeds the number of digits requested.")
		}
		
	num.zeros <- digits - num.digs
	zero.text <- paste(rep(0, num.zeros), collapse = "")
	final.num <- paste(zero.text, x, sep = "")
	return(final.num)
	
}