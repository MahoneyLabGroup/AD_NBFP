#This function finds the area under an ROC curve
smooth.auc <- function(roc.curve, plot.results = TRUE){

	require("DescTools")

	best.fit <- try(smooth.spline(roc.curve[,1], roc.curve[,2]), silent = TRUE)
	if(class(best.fit) == "try-error"){
		roc.auc <- AUC(roc.curve[,1], roc.curve[,2])
		}else{
		roc.auc <- AUC(best.fit$x, best.fit$y)
		}
	
	if(plot.results){
		plot(roc.curve[,1], roc.curve[,2], pch = 16, col = "darkgray", xlab = "False Positive Rate", ylab = "True Positive Rate")
		if(class(best.fit) != "try-error"){
			points(best.fit$x, best.fit$y, type = "l", col = "#1f78b4", lwd = 3)
			}
		abline(0,1, col = "darkgray")
		text(x = 0.7, y = 0.25, labels = paste("AUC = ", signif(roc.auc, 2)), cex = 1.5)
		}
	
	
	invisible(roc.auc)
	
}
