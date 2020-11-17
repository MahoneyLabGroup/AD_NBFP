#This function gets the variance explained from 
#a linear model

var.exp <- function(model, position){
	factor.af <- anova(model)
	factor.afss <- factor.af$"Sum Sq"
	factor.PctExp = (factor.afss/sum(factor.afss))*100
	pct.exp <- factor.PctExp[position]
	return(pct.exp)
	}