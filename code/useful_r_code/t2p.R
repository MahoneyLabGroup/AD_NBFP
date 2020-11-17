#This function generates a p value given a t value and degrees of freedom
#t.val <- seq(-4, 4, 0.1)
# df <- 1
#df <- nrow(cross$pheno)-1
#t.val <- effects1[,1]

t2p <- function(t.val, df){
	pvals <- 2*pt(-abs(t.val),df=df)
	return(pvals)
	}