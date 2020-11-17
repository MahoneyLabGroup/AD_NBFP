#This function gets the value nearest the
#given percentile in a distribution


	what.percentile <- function(distX, X){
		ext <- which(distX >= X)
		perc <- (1 - (length(ext)/length(distX)))*100
		return(perc)
		}
