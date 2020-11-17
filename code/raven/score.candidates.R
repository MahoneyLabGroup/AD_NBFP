#This function classifies positional (or other) candidates
#based on trained SVM models

score.candidates <- function(results.dir = ".", candidate.genes, verbose = FALSE){
	
	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table

	
	#go through each module
	for(i in 1:length(module.dir)){
		
		svm.csv.file <- file.path(module.dir[i], "Candidate.Gene.SVM.Scores.csv")
		fp.csv.file <- paste0(module.dir[i], "/Candidate.Gene.FP.Rates.csv")
		svm.jpg.file <- paste0(module.dir[i], "/Candidate.Gene.SVM.Scores.jpg")
		fp.jpg.file <- paste0(module.dir[i], "/Candidate.Gene.FP.Rates.jpg")
		
		if(file.exists(svm.csv.file)){next()}
		
		if(verbose){cat(dir.table[i,1], dir.table[i,2], "\n")}
		#determine the number of models that were run
		model.files <- list.files(module.dir[i], pattern = "SVM_Results")
		n.models <- length(model.files)
		
		#read in the adjacency matrix for the gene set.
		#this will be either a matrix decomposed by SVD
		#or one that has not been decomposed.
		decomp.mat.file <- list.files(module.dir[i], pattern = "decomp", full.names = TRUE)
		if(length(decomp.mat.file) == 0){
			stop(paste0("I could not find a non.decomp.mat.RData a decomp.mat.RData file in ", dir.table[i,1], " ", dir.table[i,2], ". Please make sure generate.triage.models() was run."))
			}
		decomp.mat <- readRDS(decomp.mat.file)
		

		#find the candidate genes in the adjacency matrix
		#This gives us the connections of each of the candidate genes to 
		#our true positive genes. Each SVM model uses these weights
		#to give each candidate gene a decision value
		region.gene.locale <- match(candidate.genes, colnames(decomp.mat))
		region.gene.locale <- region.gene.locale[which(!is.na(region.gene.locale))]
		test.mat <- t(decomp.mat[,region.gene.locale,drop=FALSE])

		#initialize a matrix to hold the prediction values for each gene
		#by each SVM model
		locus.gene.svm.predictions <- locus.gene.fp <- matrix(NA, ncol = nrow(test.mat), nrow = n.models)
		rownames(locus.gene.svm.predictions) <- rownames(locus.gene.fp) <- paste0("model", 1:n.models)
		colnames(locus.gene.svm.predictions) <- colnames(locus.gene.fp) <- rownames(test.mat)
		locus.gene.fp.tp <- vector(mode = "list", length = n.models)
		
		#calculate decision values for each gene in the locus
		#also calculate a false positive rate for each gene in the locus
		#this value is based on the decision values of the annotated positives
		for(md in 1:n.models){
			if(verbose){report.progress(md, n.models)}
			
			model <- readRDS(paste0(module.dir[i], "/SVM_Results", md, ".RData"))
			prediction <- predict(model$opt.model, newdata = test.mat, decision.values = TRUE)
			locus.gene.svm.predictions[md,] <- attributes(prediction)$decision.values
			locus.gene.fp.tp[[md]] <- svm.prediction.fp(attributes(prediction)$decision.values, model$opt.model)
			locus.gene.fp[md,] <- locus.gene.fp.tp[[md]][,1]
			}
		
		model.means <- colMeans(locus.gene.svm.predictions, na.rm = TRUE)
		mean.order <- order(model.means, decreasing = TRUE)
		
		fig.width = max(7, (length(candidate.genes)/33))
		fig.height = max(5, (length(candidate.genes)/66))
		
		jpeg(svm.jpg.file, width = fig.width, height = fig.height, units = "in", res = 300)
		boxplot(locus.gene.svm.predictions[,mean.order], main = "Candidate Gene SVM Scores", las = 2, cex.axis = 0.2)
		abline(h = 0)
		dev.off()
		
		jpeg(fp.jpg.file, width = fig.width, height = fig.height, units = "in", res = 300)
		fp.means <- colMeans(locus.gene.fp, na.rm = TRUE)
		mean.order <- order(fp.means, decreasing = FALSE)
		boxplot(locus.gene.fp[,mean.order], main = "Candidate Gene FP Rates", las = 2, cex.axis = 1)
		abline(h = 0)
		dev.off()
		
		write.table(locus.gene.svm.predictions, svm.csv.file, sep = ",", quote = FALSE)
		
		write.table(locus.gene.fp, fp.csv.file, sep = ",", quote = FALSE)
		
		}
	
	
}
