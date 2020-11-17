iterate.function <-
function(dir.list = "dir.list.txt", FUN){

	base.dir <- getwd()
	
	print(paste("working in base directory:", base.dir))

	experiment.list = read.table(dir.list, stringsAsFactors = FALSE, header = TRUE)
	
	exp.dir.names = experiment.list[,1]

	result <- list()
	no.results.exp <- NULL
	result.counter <- 1
	for(e in 1:length(exp.dir.names)){
		print(paste(exp.dir.names[e]), sep = "")
		setwd(paste(base.dir, exp.dir.names[e], sep = "/"))
			FUN <- match.fun(FUN)
			results <- FUN()
			if(length(results) > 0){
				result[[result.counter]] <- results
				result.counter <- result.counter + 1
				}else{
					no.results.exp <- c(no.results.exp, e)
					}
		}
		
		if(length(no.results.exp) > 0){
			exp.dir.names <- exp.dir.names[-no.results.exp]
			}
		if(length(result) > 0){
			names(result) <- exp.dir.names
			return(result)
			}

	}
