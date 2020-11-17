#This function draws a square
#on an existing plot between
#a minimum xy coordinate and
#a maximum xy coordinate



draw.rectangle <- function(min.x, max.x, min.y, max.y, border.col = "black", fill = NULL, lwd = 1){
	
	# plot.new()
	# plot.window(c(1,10), c(1,10))
	# min.x = 1; max.x = 5; min.y = 7; max.y = 9
	polygon(x = c(min.x, max.x, max.x, min.x), y = c(max.y, max.y, min.y, min.y), border = border.col, col = fill, lwd = lwd)
	 
	
}