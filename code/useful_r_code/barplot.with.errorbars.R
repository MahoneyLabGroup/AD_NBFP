barplot.with.errorbars <- function(data.matrix, margin = 2, error.bars = c("sd", "se"), ...){
	
	if(margin == 1){
		data.matrix <- t(data.matrix)	
		}	
	
	plotting.args <- list(...)
	# print(plotting.args)

	bar.heights <- colMeans(data.matrix)

	if(error.bars == "se"){
		bar.error <- bar.error/sqrt(dim(data.matrix)[1])
		}else{
		bar.error <- apply(data.matrix, 2, sd)			
			}
	
	
	y.min <- min(bar.heights-bar.error); y.max <- max(bar.heights+bar.error)
	plotting.args$height <- bar.heights
	plotting.args$ylim <- c(y.min, y.max)
	#get the midpoints of the bars
	x.pos <- barplot(bar.heights, ylim = c(y.min, y.max), plot = FALSE)
	
	
	#draw the plot with the list of arguments
	do.call(barplot, plotting.args)
	x.range <- max(x.pos) - min(x.pos)
	segments(x0 = x.pos[,1], y0 = bar.heights-bar.error, x1 = x.pos[,1], y1 = bar.heights+bar.error) #the vertical bars
	segments(x0 = (x.pos[,1] - (x.range*0.01)), y0 = bar.heights+bar.error, x1 = (x.pos[,1] + (x.range*0.01)), y1 = bar.heights+bar.error) #top horizontal bars
	segments(x0 = (x.pos[,1] - (x.range*0.01)), y0 = bar.heights-bar.error, x1 = (x.pos[,1] + (x.range*0.01)), y1 = bar.heights-bar.error) #bottom horizontal bars
	
	
}
