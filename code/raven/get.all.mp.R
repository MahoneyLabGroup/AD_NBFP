#This function gets the list of all mp terms from mice

get.all.mp <- function(){

	require(InterMineR)
	require(biomaRt)

	all.var <- ls(globalenv())
	
	mine <- initInterMine("http://www.mousemine.org/mousemine")
	# getTemplates(mine)
	mp.query <- getTemplateQuery(mine, "Lookup_MPhenotype")
	mp.query$where[[2]]$value <- "*"

	mp.terms <- runQuery(mine, mp.query)	
	
	mp.genes <- vector(mode = "list", length = nrow(mp.terms))
	names(mp.genes)  <- mp.terms[,1]

	for(i in 1:nrow(mp.terms)){	
		print(i)
		mp.genes[[i]] <- get.mp.genes(mp.terms[i,1], "mouse", verbose = FALSE)
		}


	return(mp.genes)
	
	}