#This function removes files with a given pattern from 
#the module directories in a TRiAGE project

remove.file <- function(results.dir, pattern){
	module.dirs <- get.module.dir(results.dir)
	for(i in 1:length(module.dirs)){
		to.remove <- list.files(path = module.dir.info[i], pattern = pattern, full.names = TRUE)
		if(length(to.remove) > 0){
			for(j in 1:length(to.remove)){
				system(paste("rm", to.remove[j]))
				}
			}
		}
	}