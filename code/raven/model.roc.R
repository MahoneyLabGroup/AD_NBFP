#calcuate ROC curves for all models	in
#a directory.
#these tell us how separable our term
#genes are from a balanced set of randomly
#selected genes outside the term. 
#positives are annotated positives
#and negatives are the randomly selected
#genes outside the term used in the models

model.roc <- function(path = "."){
	
	all.files <- list.files(path = path, pattern = "SVM_Results")
	
	all.roc <- lapply(1:trait.genes.svm.trials, function(x) svm.model.roc(readRDS(paste0("SVM_Results", x, ".RData"))$opt.model))

	pdf("ROC.Model.pdf")
	model.auc <- smooth.auc(roc.curve = Reduce("rbind", all.roc), plot.results = TRUE)
	dev.off()
	}