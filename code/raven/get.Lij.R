#this function returns the cape locus network
#this is the marker-marker network without covariates

get.Lij <- function(data.obj, collapsed.net = TRUE, remove.covar = TRUE){
		
	if(collapsed.net){
		net <- data.obj$collapsed.net
		}else{
		net <- data.obj$full.net	
		}
	
	covar.info <- get.covar(data.obj)
	covar.locale <- which(rownames(net) %in% covar.info$covar.names)
	
	just.markers <- net[,1:nrow(net)]
	if(remove.covar && length(covar.locale) > 0){
		final.net <- just.markers[-covar.locale,-covar.locale]
		}else{
		final.net <- just.markers	
		}
	
	return(final.net)
}