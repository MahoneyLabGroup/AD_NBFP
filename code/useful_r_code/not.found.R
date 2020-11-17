#This function takes in two vectors
#and reports which items in V1 cannot
#be found in V2

not.found <- function(V1, V2){
	
	only.in1 <- setdiff(V1, V2)
	cat(only.in1, sep = "\n")	
	
}