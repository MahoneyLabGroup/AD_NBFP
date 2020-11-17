delete.files <-
function(dir.list, pattern){

	curr.dir <- getwd()
	dirs <- read.table(dir.list, stringsAsFactors = FALSE, header = TRUE)	
	
	for(i in 1:length(dirs[,1])){
		setwd(dirs[i,1])
		unlink(paste("*", pattern, "*", sep = ""))
		setwd(curr.dir)	
		}

	}
