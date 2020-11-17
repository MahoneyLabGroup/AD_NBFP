#This function gets the value nearest the
#given percentile in a distribution


	get.percentile <- function(distX, percentile){
		ranked <- rank(distX)
		#find the position that the percentile should be in the ranked list
		perc.pos <- round(length(distX)*(percentile/100))
		#find the closest thing to this position 
		locale <- which(abs((ranked - perc.pos)) == min(abs(ranked - perc.pos)))[1]
		return(distX[locale])
		}
