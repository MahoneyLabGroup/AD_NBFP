#This function performs a mediation analysis.
#It asks for the variance of the outcome variable
#(out.var) explained by the explanatory variable 
#(exp.var) using a linear model (lm). It then places 
#each of the mediators in the linear model, and 
#gets the new variance explained by the exp.var
#after mediating on each mediator.
#out.var is a matrix of n samples in rows and
#v columns corresponding to different variables,
#for example fear conditioning, amyloid, etc.
#exp.var is a vector of length n samples. 
#mediators is a matrix with n rows and 
#m mediators
#This function does not take covariates into account
#so all variables must be adjusted for all covariates
#before applying this function.

#library(rgl)
#exp.var <- rnorm(100)
#out.var <- matrix(rnorm(300, sd = 0.7), nrow = 100, ncol = 3)
#out.var <- apply(out.var, 2, function(x) (exp.var*runif(1,0,2))+x)
#colnames(out.var) <- paste0("trait", 1:ncol(out.var))
#mediators <- matrix(rnorm(10000), nrow = 100, ncol = 100)
# mediators <- apply(mediators, 2, function(x) (exp.var*out.var[,runif(1, 1, 3)]-x*0.5))
# colnames(mediators) <- paste0("mediator", 1:ncol(mediators))
#plot3d(exp.var, out.var[,1], mediators[,1])
#if calc.mediator.exp.var is TRUE, this function will also calculate the
#percent variance explained by each mediator alone without the explanatory variable.

mediation.analysis <- function(out.var, exp.var, mediators, plot.results = FALSE, 
calc.mediator.exp.var = FALSE){
	
	#========================================================
	#do some quick checks for variance in the input matrices
	#========================================================
	
	test.var <- round(apply(out.var, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		stop("some outcome variables have 0 variance.")
		}
	
	test.var <- round(apply(mediators, 2, function(x) var(x, na.rm = TRUE)), 2)
	if(any(test.var == 0)){
		zero.locale <- which(test.var == 0)
		mediators <- mediators[,-zero.locale]
		warning(paste("removing", length(zero.locale), "mediators with zero variance"))
		}
	
	test.var <- round(apply(mediators, 1, function(x) var(x, na.rm = TRUE)), 4)
	na.locale <- which(is.na(test.var))
	zero.locale <- which(test.var == 0)
	if(length(na.locale) > 0 || length(zero.locale) > 0){
		remove.locale <- unique(c(na.locale, zero.locale))
		mediators <- mediators[-remove.locale,]
		out.var <- out.var[-remove.locale,]
		exp.var <- exp.var[-remove.locale]
		warning(paste("removing", length(remove.locale), "individuals with zero expression variance"))
		}

	
	#========================================================
	#find the variance of the outcome measure explained by each
	#explanatory variable
	#========================================================
		
	# boxplot(out.var[,2]~exp.var)
	factor.test <- apply(out.var, 2, function(x) lm(x~exp.var))
	orig.pct.exp <- unlist(lapply(factor.test, function(x) var.exp(x, 1)))
	
	#========================================================
	#put each mediator in the model and calculate the new
	#variance explained
	#========================================================	
	all.var.exp <- just.var.exp <- matrix(NA, nrow = ncol(mediators), ncol = ncol(out.var))
	colnames(all.var.exp) <- colnames(just.var.exp) <- colnames(out.var)
	rownames(all.var.exp) <- rownames(just.var.exp) <- colnames(mediators)

	
	for(i in 1:ncol(out.var)){
		cat("calculating mediators of", colnames(out.var)[i], "...\n")
		mediator.models <- apply(mediators, 2, function(x) lm(out.var[,i]~x+exp.var))
		all.var.exp[,i] <- unlist(lapply(mediator.models, function(x) var.exp(x, 2)))

		if(calc.mediator.exp.var){
			just.var <- apply(mediators, 2, function(x) lm(out.var[,i]~x))
			just.var.exp[,i] <- unlist(lapply(just.var, function(x) var.exp(x, 1)))
		}
	}
	
	
	if(plot.results){
		xlim = c(0, (ncol(out.var)+1))
		ylim = c(0, ceiling(max(c(all.var.exp, orig.pct.exp))))
		plot.new()
		plot.window(xlim = xlim, ylim = ylim)

		#add lines for the original variance explained by the explanatory 
		#variable for each outcome variable
		segments(x0 = 1:ncol(out.var)-0.2, y0 = orig.pct.exp, x1 = 1:ncol(out.var)+0.2, lwd = 3, col = "#3182bd")

		#add boxplots for the new variance explained after the mediators have been added
		for(i in 1:ncol(all.var.exp)){
			boxplot(all.var.exp[,i], at = i, add = TRUE, axes = FALSE)
			}
		axis(2)		
		mtext("Percent Variance Explained", side = 2, line = 2.5)
		par(xpd = TRUE)
		text(x = 1:ncol(out.var), y = 0-(ylim[2]*0.05), labels = colnames(out.var))	
		par(xpd = FALSE)
	}
	
	if(calc.mediator.exp.var){
		final.result <- list("mediation.var.exp" = rbind(orig.pct.exp, all.var.exp), 
		mediator.var.exp = just.var.exp)
		}else{
		final.result <- rbind(orig.pct.exp, all.var.exp)
		}
	
	return(final.result)
	
	}