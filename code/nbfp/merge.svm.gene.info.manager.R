#This function 
merge.svm.gene.info.manager <- function(results.dir = ".", gene.info.table){
	
	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table

	all.results <- vector(mode = "list", length = length(module.dir))
	names(all.results) <- dir.table[,2]
	
	for(i in 1:length(module.dir)){
		
		svm.csv.file <- paste0(module.dir[i], "/Candidate.Gene.SVM.Scores.csv")
		fp.csv.file <- paste0(module.dir[i], "/Candidate.Gene.FP.Rates.csv")
		results.file <- paste0(module.dir[i], "/Candidate.Gene.Results.csv")
		plot.file <- paste0(module.dir[i], "/Candidate.Gene.SVM.Results.jpg")
		
		svm.scores <- read.csv(svm.csv.file, stringsAsFactors = FALSE)
		mean.svm <- colMeans(svm.scores)
		gene.ids <- gsub("X", "", colnames(svm.scores))
		
		fp.rates <- read.csv(fp.csv.file, stringsAsFactors = FALSE)
		mean.fp <- colMeans(fp.rates)

		common.genes <- intersect(as.character(t(gene.info.table[,"GENE"])), gene.ids)
		gene.locale.table <- match(common.genes, as.character(t(gene.info.table[,"GENE"])))
		# head(cbind(gene.ids, gene.info.table[gene.locale.table,]))
		gene.locale.svm <- match(common.genes, gene.ids)
		# head(cbind(common.genes, gene.ids[gene.locale.svm]))
		
		final.table <- cbind(gene.info.table[gene.locale.table,], mean.svm[gene.locale.svm], mean.fp[gene.locale.svm])
		colnames(final.table)[c(10,11)] <- c("Mean.SVM.Score", "Mean.FP.Rate") #Should not be hard-coded since columns may change depending on input table
		
		write.table(final.table, results.file, sep = ",", quote = FALSE, row.names = FALSE)
		
		mean.gene.position <- rowMeans(final.table[,c("START", "STOP")])
		
		jpeg(plot.file, height = 7, width = 10, units = "in", res = 300)
		plot.new()
		plot.window(xlim = c(min(mean.gene.position), max(mean.gene.position)), ylim = c(min(final.table[,"Mean.SVM.Score"]), max(final.table[,"Mean.SVM.Score"])))
		text(x = mean.gene.position, y = final.table[,"Mean.SVM.Score"], labels = final.table[,"GENE"], cex = 0.7)
		axis(1);axis(2)
		mtext("Genomic Position", side = 1, line = 2.5)				
		mtext("SVM Score", side = 2, line = 2.5)
		abline(h = 0)
		dev.off()
		
		all.results[[i]] <- final.table
		}
	invisible(all.results)

}