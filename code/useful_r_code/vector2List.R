#This function turns a vector V into a list
#based on categories in the vector cats

vector2List <- function(V, cats){
	u_cats <- sort(unique(cats))
	vlist <- lapply(u_cats, function(x) V[which(cats == x)])
	names(vlist) <- u_cats
	return(vlist)
	}