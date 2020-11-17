#This is a general version of the order.strains() function 
#used in the Epigenetics projects.
#given two vectors of names, it returns
#the order the second vector should be put in
#to match the first. It expects key, which
#is stored in the file strain.color.table.txt
#If this file is changed, this function will need
#to change

match.order <- function(Vfixed, Vreorder, key){
	
	#figure out which column each strain vector matches the best
	match.col <- function(strainV, key){
		strain.matches <- apply(key, 2, function(x) match(strainV, x))
		num.matches <- apply(strain.matches, 2, function(x) length(which(!is.na(x))))
		best.match <- which.max(num.matches)
		return(as.vector(best.match))
		}
	
	Vfixed.col <- match.col(Vfixed, key)
	Vreorder.col <- match.col(Vreorder, key)
	
	#sort key to match the order of the Vfixed vector
	fixed.order <- match(Vfixed, key[,Vfixed.col])
	#cbind(Vreorder, key[fixed.order,])
	fixed.order.table <- key[fixed.order,]

	#Then reorder the Vreorder vector to match that vector
	reorder.order <- match(fixed.order.table[,Vreorder.col], Vreorder)
	#cbind(Vreorder[reorder.order], Vfixed)

	return(reorder.order)

}