#This function sources functions in listed
#directories

source.fun <- function(fun.dir){
	
	for(f in 1:length(fun.dir)){
		all.fun <- list.files(path = fun.dir[f], pattern = ".R", full.names = TRUE)
		for(a in 1:length(all.fun)){source(all.fun[a])}
		}

}