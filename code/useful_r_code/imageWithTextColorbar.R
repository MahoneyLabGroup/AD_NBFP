#testing: 
# mat <- matrix((1:100), 10, 10)
# mat <- matrix(rnorm(100), 10, 10)
#orientation = c("h", "v"); num.labels = 10; cex = 0.5; split.at.vals = FALSE; split.points = 0
#col.scale = c("green", "purple", "orange", "blue", "gray"); light.dark = "f"; class.mat = NULL
#grad.dir = c("high", "low", "middle", "ends"); color.fun = c("linear", "exponential")
#exp.steepness = 1; global.color.scale = FALSE; global.min = NULL; global.max = NULL 
#axis.line = -3; use.pheatmap.colors = TRUE

imageWithTextColorbar <- function(mat, orientation = c("h", "v"), num.labels = 10, cex = 0.5, split.at.vals = FALSE, split.points = 0, col.scale = c("green", "purple", "orange", "blue", "gray"), light.dark = "f", class.mat = NULL, grad.dir = c("high", "low", "middle", "ends"), color.fun = c("linear", "exponential"), exp.steepness = 1, global.color.scale = FALSE, global.min = NULL, global.max = NULL, axis.line = -3, use.pheatmap.colors = FALSE, ax.min = NULL, ax.max = NULL, n.ax.ticks = NULL,
hadj = NA, padj = NA, bounding.box = TRUE){

	
		library(grid)
	
		if(length(light.dark) < length(col.scale)){light.dark <- rep(light.dark, length(col.scale))}
		
		orientation.check <- grep("v", orientation)
		if(length(orientation.check) > 0){
			orientation <- "v"
			}
		
		get.default.col.fun <- grep("lin", color.fun)
		if(length(get.default.col.fun) > 0){
			color.fun = "linear"
			}
		
		end.fudge.factor = 10^-10

		if(is.null(class.mat)){
			class.mat <- matrix(1, dim(mat)[1], dim(mat)[2])
			}

		if(split.at.vals){
			for(p in 1:length(split.points)){
				class.mat[which(mat >= split.points[p])] <- p+1
				}
			class.boundaries <- c(min(mat, na.rm = TRUE), split.points, max(mat, na.rm = TRUE))
			}else{
			split.points <- NULL	
			class.boundaries <- c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE))
			}
			
		

		num.classes <- length(unique(as.vector(class.mat[which(!is.na(class.mat))])))
		if(num.classes == 1){
			class.mat <- matrix(1, dim(mat)[1], dim(mat)[2])
			}
			
		class.cols <- col.scale[1:num.classes]

			
		get.default <- grep("h", grad.dir)
		if(length(get.default) > 0){
			grad.dir <- "high"
			}
		
		dir.list <- vector(mode = "list", length = num.classes)
		if(grad.dir == "high"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- 1:3
				}
			}
		if(grad.dir == "low"){
			for(i in 1:length(dir.list)){
				dir.list[[i]] <- 3:1
				}
			}
		if(grad.dir == "middle"){
			if(length(dir.list) != 2){stop("I can only color the middle if there are exactly two classes")}
			dir.list[[1]] <- 1:3
			dir.list[[2]] <- 3:1
			}
			
			
		if(grad.dir == "ends"){
			if(length(dir.list) != 2){stop("I can only color the ends if there are exactly two classes")}
			dir.list[[1]] <- 3:1
			dir.list[[2]] <- 1:3
			}
	
	
		num.classes = length(unique(as.vector(class.mat)))

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
		fill.color.ramp <- function(mat, class.mat, global){
			full.color.key <- NULL
			num.col <- NULL
			
			#make color scales for each class or globally as defined
			color.scales <- vector(mode = "list", length = num.classes)

			ColorRamp <- matrix(NA, dim(mat)[1], dim(mat)[2])
			
			for(cl in 1:num.classes){
				if(global){
					if(is.null(global.min)){min.cl = min(mat, na.rm = TRUE)}else{min.cl = global.min}
					if(is.null(global.max)){max.cl = max(mat, na.rm = TRUE)}else{max.cl = global.max}	
					}else{
					min.cl <- class.boundaries[cl]
					max.cl <- class.boundaries[cl+1]
					}
					
				if(color.fun == "linear"){
					ColorLevels <- seq(min.cl, max.cl, length=256)
					}else{
					ColorLevels <- exp.color.fun(min.cl, max.cl, steepness = exp.steepness, num.cols=256)	
					}

				#make the function to generate 
				col.vals <- get.color(col.scale[cl], light.dark[cl])
				color.scales[[cl]] <- colorRampPalette(col.vals[dir.list[[cl]]])
				
				#find the entries in each class
				entry.locale <- which(class.mat == cl)
				entry.vals <- mat[entry.locale]
				color.key <- cbind(ColorLevels, do.call(color.scales[[cl]], list(256)))
				entry.cols <- bin.cols(color.key, entry.vals)
				entry.order <- order(entry.locale)
				ColorRamp[entry.locale[entry.order]] <- entry.cols
				if(split.at.vals){
					if(cl == 1){
						split.locale <- which(as.numeric(color.key[,1]) <= split.points[cl])
						}
					if(cl > 1){
						split.locale <- intersect(which(as.numeric(color.key[,1]) <= split.points[cl]), which(as.numeric(color.key[,1]) > split.points[cl-1]))
						}
					if(length(split.locale) > 0){
						full.color.key <- rbind(full.color.key, color.key[split.locale,])
						num.col <- c(num.col, dim(color.key)[1])
						}else{
						split.locale <- which(as.numeric(color.key[,1]) > split.points[cl-1])
						if(length(split.locale) > 0){
							full.color.key <- rbind(full.color.key, color.key[split.locale,])
							num.col <- c(num.col, dim(color.key)[1])
							}
						}
					}else{
					full.color.key <- color.key	
					num.col <- dim(color.key)[1]
					}
				}
			
			return(list(full.color.key, num.col))
			
			}
		
		#============================================================================


		if(use.pheatmap.colors){
			pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
			bks <- segment.region(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), 
			length(pal), alignment = "ends")
			col.key <- cbind(bks[1:length(pal)], pal)
			num.cols <- nrow(col.key)
			}else{	
			col.key.list <- fill.color.ramp(mat, class.mat, global.color.scale)
			col.key <- col.key.list[[1]]
			num.cols <- col.key.list[[2]]
			}

		# print(num.cols)
		# print(dim(col.key))

	# par(mfrow = c(1,2))
		# plot(x = as.numeric(col.key[,1]), y = rep(1, dim(col.key)[1]), col = col.key[,2], pch = "|", cex = 1, xlab = "", ylab = "")
		# abline(v = split.points)

		if(!global.color.scale){
			plot.new()
			if(orientation == "h"){
				plot.window(xlim = c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), ylim = c(0,1))
				axis.side = 1
				}else{
				plot.window(ylim = c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), xlim = c(0,1))	
				axis.side = 2
				}
			par(mar = c(3,0,0,0))
			start.x <- min(mat, na.rm = TRUE)
			# for(i in 2:100){
			for(i in 2:dim(col.key)[1]){
				if(orientation == "h"){
				polygon(x = c(start.x, rep(as.numeric(col.key[i,1]), 2), start.x), y = c(0,0,1,1), col = col.key[i,2], border = NA)
				}else{
				polygon(y = c(start.x, rep(as.numeric(col.key[i,1]), 2), start.x), x = c(0,0,1,1), col = col.key[i,2], border = NA)
				}
				start.x <- as.numeric(col.key[i, 1])
				}

			if(!is.null(n.ax.ticks)){ #if a number of tick marks is specified
				if(is.null(ax.min)){ax.min <- floor(min(mat, na.rm = TRUE))}
				if(is.null(ax.max)){ax.max <- ceiling(max(mat, na.rm = TRUE))}
				at <- segment.region(ax.min, ax.max, n.ax.ticks, alignment = "ends")
				}else{
				at = NULL #otherwise use the default axis
				}
			axis(axis.side, line = axis.line, cex.axis = cex, at = at, padj = padj, hadj = hadj)

			if(bounding.box){
				if(orientation == "h"){
					polygon(x = c(rep(min(mat, na.rm = TRUE), 2), rep(max(mat, na.rm = TRUE),2)), y = c(0,1,1,0))
					}else{
					polygon(y = c(rep(min(mat, na.rm = TRUE), 2), rep(max(mat, na.rm = TRUE),2)), x = c(0,1,1,0))	
					}
				}
			}else{
			num.mat <- as.matrix(as.numeric(col.key[,1]), ncol = 1)
			if(orientation == "h"){
				image(x = as.vector(num.mat), y = 1, z = num.mat, col= col.key[,2], xlab="",ylab="", yaxt = "n", cex.axis = cex)
				}else{
				# Original
				image(x = 1, y = as.vector(num.mat), z = t(num.mat), col= col.key[,2], xlab="",ylab="",xaxt="n", cex.axis = cex)
				}
			}
					


	}