distribute.file <-
function(dir.list, base.name){

	base.dir <- get.base.dir()
	setwd(base.dir)
	dirs <- read.table(dir.list, stringsAsFactors = FALSE, header = TRUE)	
	
	
	for(i in 1:length(dirs[,1])){
		system(paste("cp", base.name, dirs[i,"dir.label"]))
		}
	
	}
