#This function is based on imageWithText
#it gets colors for numeric values

colors.from.values <- function(vals, split.at.vals = FALSE, split.points = 0, 
col.scale = c("green", "purple", "orange", "blue", "brown", "gray"), light.dark = "f", 
grad.dir = c("high", "low", "middle", "ends"), color.fun = c("linear", "exponential"), 
exp.steepness = 1, global.color.scale = FALSE, global.min = NULL, global.max = NULL, 
use.pheatmap.colors = FALSE, na.col = "lightgray"){

		require(grid)
	 	class.mat = NULL
		#make sure Inf and -Inf are coded as NA
		vals[which(!is.finite(vals))] <- NA
		
		if(length(which(is.na(vals))) == length(vals)){
			return()
			}

		
		get.default.col.fun <- grep("lin", color.fun)
		if(length(get.default.col.fun) > 0){
			color.fun = "linear"
			}
		
		end.fudge.factor = 10^-10

		if(is.null(class.mat)){
			class.mat <- rep(1, length(vals))
			}

		if(split.at.vals){
			for(p in 1:length(split.points)){
				class.mat[which(vals >= split.points[p])] <- class.mat[which(vals >= split.points[p])] + 1
				}
			}else{
			split.points <- NULL	
			}
		class.mat[which(is.na(vals))] <- NA

		while(min(class.mat, na.rm = TRUE) > 1){class.mat <- class.mat - 1}
		
		classes <- sort(unique(as.vector(class.mat)))
		num.classes <- length(classes)
		if(num.classes == 1){
			class.mat <- rep(1, length(vals))
			}
	
		if(length(col.scale) == (length(split.points)+1)){
			class.cols <- col.scale
			}else{
			class.cols <- col.scale[classes]
			}
		if(length(col.scale) < num.classes){
			extra.cols <- num.classes - length(col.scale)
			col.scale <- c(col.scale, rep(col.scale, ceiling(extra.cols/length(col.scale))))
			}
			
		get.default <- grep("h", grad.dir)
		if(length(get.default) > 0){
			grad.dir <- "high"
			}
		
		max.col = 4
		dir.list <- vector(mode = "list", length = num.classes)
		names(dir.list) <- classes
		if(grad.dir == "high"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- 1:max.col
				}
			}
		if(grad.dir == "low"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- max.col:1
				}
			}
		if(grad.dir == "middle"){
			if(length(dir.list) != 2){stop("I can only color the middle if there are exactly two classes")}
			dir.list[[1]] <- 1:max.col
			dir.list[[2]] <- max.col:1
			}
			
			
		if(grad.dir == "ends"){
			if(length(dir.list) != 2){stop("I can only color the ends if there are exactly two classes")}
			dir.list[[1]] <- max.col:1
			dir.list[[2]] <- 1:max.col
			}



		#============================================================================
		#internal functions
		#============================================================================
		#This function takes in a matrix of values matched with colors, and 
		#a vector of values. It matched up the appropriate color for each
		#value in the vector
		bin.cols <- function(color.key, V){
			color.v <- rep(NA, length(V))
			for(i in 1:length(V)){
				diff.v  <- V[i] - as.numeric(color.key[,1])
				closest.val <- which(abs(diff.v) == min(abs(diff.v)))[1]
				color.v[i] <- color.key[closest.val,2]
				}
			return(color.v)
			}

		#This function generates the matrix of colors to use in 
		#the raster function
		fill.color.ramp <- function(vals, class.mat, global){

			#make color scales for each class or globally as defined
			color.scales <- vector(mode = "list", length = num.classes)
			names(color.scales) <- classes

			ColorRamp <- rep(NA, length(vals))
			num.classes = length(unique(as.vector(class.mat)))
			
			for(cl in 1:length(classes)){
				if(global.color.scale){
					if(is.null(global.min)){min.cl = min(vals, na.rm = TRUE)}else{min.cl = global.min}
					if(is.null(global.max)){max.cl = max(vals, na.rm = TRUE)}else{max.cl = global.max}	
					}else{
					if(length(which(class.mat == classes[cl])) > 0 && !all(is.na(vals[which(class.mat == classes[cl])]))){
						min.cl <- min(vals[which(class.mat == classes[cl])], na.rm = TRUE)
						max.cl <- max(vals[which(class.mat == classes[cl])], na.rm = TRUE)
						}else{
						min.cl <- NA	
						}
						
					}
					
				if(!is.na(min.cl)){
					if(color.fun == "linear"){
						ColorLevels <- seq(min.cl, max.cl, length=256)
						}else{
						ColorLevels <- exp.color.fun(min.cl, max.cl, steepness = exp.steepness, num.cols=256)	
						}
	
					#make the function to generate 
					col.vals <- get.color(col.scale[cl], light.dark)
					color.locale <- which(names(color.scales) == classes[cl])
					color.scales[[color.locale]] <- colorRampPalette(col.vals[dir.list[[color.locale]]])
					
					#find the entries in each class
					entry.locale <- which(class.mat == classes[cl])
					if(length(entry.locale) > 0){
						entry.vals <- vals[entry.locale]
						color.key <- cbind(ColorLevels, do.call(color.scales[[color.locale]], list(256)))
						entry.cols <- bin.cols(color.key, entry.vals)
						entry.order <- order(entry.locale)
						ColorRamp[entry.locale[entry.order]] <- entry.cols
						}
					}
				}
			
			return(ColorRamp)
			
			}
		
		#============================================================================



		if(use.pheatmap.colors){
			pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
			bks <- pheatmap:::generate_breaks(vals, length(pal), center = F)
			ColorRamp <- pheatmap:::scale_colours(vals, col=pal, breaks=bks, na_col = na.col)
		}else{
			ColorRamp = fill.color.ramp(vals, class.mat, global.color.scale)
		}
		 	
		zmin <- min(vals, na.rm = TRUE); zmax <- max(vals, na.rm = TRUE)

		na.locale <- which(is.na(vals))
		if(length(na.locale) > 0){
			vals[na.locale] <- 0
			}
		
		all.ind <- which(!is.na(vals), arr.ind = TRUE)
		if(length(na.locale) > 0){
			vals[na.locale] <- NA
			}
		
		#translate the ColorRamp matrix to rgb matrices so we can use grid.
		rgb.mat <- col2rgb(ColorRamp)
		col <- rgb(rgb.mat["red",]/256,rgb.mat["green",]/256,rgb.mat["blue",]/256)

		return(col)
		

	}