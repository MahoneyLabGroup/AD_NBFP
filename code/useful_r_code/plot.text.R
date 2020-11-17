plot.text <- function(the.message, srt = 0, cex = 1, x = 0.5, y = 0.5){
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	text(x, y, the.message, srt = srt, cex = cex)
	}