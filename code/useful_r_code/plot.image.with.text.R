
plot.image.with.text <- function(mat, xlab = "", ylab = "", main = NULL){

 		my.palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")
		zmin <- min(mat, na.rm = TRUE); zmax <- max(mat, na.rm = TRUE)

		na.locale <- which(is.na(mat))
		if(length(na.locale) > 0){
			mat[na.locale] <- 0
			}
		
		all.ind <- which(!is.na(mat), arr.ind = TRUE)		
		
		image(x = 1:dim(mat)[2], y = 1:dim(mat)[1], z = rotate.mat(mat), col = my.palette(50), axes = FALSE, zlim = c(zmin, zmax), main = main, xlab = xlab, ylab = ylab)

		text(all.ind[,2], rev(all.ind[,1]), labels = signif(as.vector(mat), 2), cex = 0.4)
	}