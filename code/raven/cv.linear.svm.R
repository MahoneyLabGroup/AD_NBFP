cv.linear.svm <- function(data.mat, data.labels, C.list = 4^seq(-2, 2, 1), scale = TRUE, class.weights = NULL, verbose = FALSE){
	
	require("e1071")
	
	# Initialize output data
	cv.data = NULL
	cv.data$C = C.list
	cv.data$ave.accuracy = matrix(0, nrow = length(C.list), ncol = 1)
	cv.data$opt.model = NULL
	
	# Loop and compute cross validation errors and find optimal model
	cv.data$best.cv.acc = 0
	
	for(i in 1:length(C.list)){
		
		if(verbose){
		cat("C =", C.list[i], ":")
		}
		
		# Fit model

		# curr.svm = svm(data.mat, data = data.labels, kernel = "linear", cost = C.list[i], type = "C", cross = 10, scale = scale, class.weights = class.weights)

		curr.svm = svm(x = data.mat, y = data.labels, kernel = "linear", cost = C.list[i], type = "C", cross = 10, scale = scale, class.weights = class.weights)

		
		# Compute average accuracy on folds
		curr.cv.acc = mean(curr.svm$accuracies)
		if(verbose){
		cat("\taccuracy:", curr.cv.acc, "\n")
		}
		
		cv.data$ave.accuracy[i] = curr.cv.acc
		
		# Update model
		if(curr.cv.acc > cv.data$best.cv.acc){
			cv.data$opt.model = curr.svm
			cv.data$best.cv.acc = curr.cv.acc
			}
		}
	
	# Return
	return(cv.data)
}