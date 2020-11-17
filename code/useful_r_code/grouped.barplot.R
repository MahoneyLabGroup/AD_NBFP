grouped.barplot <- function(height.list, ylim = NULL){
	
	xlim <- c(0, length(height.list))
	
	if(is.null(ylim)){ylim <- c(0, max(unlist(height.list)))	}

	
	par(xpd = TRUE, mar = c(10, 4, 4, 4))
	plot.new()
	plot.window(xlim = xlim, ylim = ylim)
	for(i in 1:length(height.list)){
		if(length(height.list[[i]]) > 1){
			x.min <- i-0.2; x.max <- i+0.2
			x.coord <- segment.region(x.min, x.max, length(height.list[[i]]), "center")			
			}else{
			x.coord <- i
			}
		points(x <- x.coord, y = sort(height.list[[i]]), lwd = 3, type = "h", col = "gray")
		text(x = i, y = 0-(ylim[2]*0.02), labels = names(height.list)[i], srt = 90, adj = 1)
		}
	axis(2)
	par(xpd = FALSE)
}