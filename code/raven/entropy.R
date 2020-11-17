#This script calcluates the entropy of a probability distribution


entropy <- function(prob.dist){
		
		#This calculation corrects for 0's
		#If there are 0's in the distribution
		#we first take the log of the distribution
		#and then replace all NaN's with 0
		zeroes <- which(prob.dist == 0)
		if(length(zeroes) > 0){
			log.dist <- log(prob.dist)
			log.dist[zeroes] <- 0
			h <- -1 * sum(prob.dist*log.dist)
			}else{
				h <- -1*sum(prob.dist*log(prob.dist))
				}
				
		return(h)

	}