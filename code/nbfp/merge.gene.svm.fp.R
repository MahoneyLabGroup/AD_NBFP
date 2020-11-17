#This function merges the gene information and SVM scores and FP Rates

merge.gene.svm.fp <- function(results.dir = ".", gene.info.table, gene.column.name = c("GENEentrezID", "GENE", "entrezID")){
  
  #directories (also works if there are multiple modules)
  module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
  module.dir <- module.dir.info$module.dir
  dir.table <- module.dir.info$dir.table
  all.results <- vector(mode = "list", length = length(module.dir))
  names(all.results) <- dir.table[,2]
  
  for(i in 1:length(module.dir)){
    svm.csv.file <- paste0(module.dir[i], "/Candidate.Gene.SVM.Scores.csv")
    fp.csv.file <- paste0(module.dir[i], "/Candidate.Gene.FP.Rates.csv")
    results.file <- paste0(module.dir[i], "/Candidate.Gene.Results.rds")
    
    #SVM scores
    svm.scores <- read.csv(svm.csv.file, stringsAsFactors = FALSE) #dataframe of 100 obs of all genes
    mean.svm <- as.data.frame(colMeans(svm.scores)) #"named num" convert to df, column names are genes, and numbers are means
    mean.svm$GENE_entrezID <- as.character(gsub("X", "", rownames(mean.svm))) #entrez gene ID column added
    colnames(mean.svm)[1] <- "mean.svm.score"
      
    #FP rates
    fp.rates <- read.csv(fp.csv.file, stringsAsFactors = FALSE)
    mean.fp <- as.data.frame(colMeans(fp.rates))
    mean.fp$GENE_entrezID <- as.character(gsub("X", "", rownames(mean.fp)))
    colnames(mean.fp)[1] <- "mean.fp.rate"
    
    #combining gene.info.table with svm scores and fp rates
    final.table <- inner_join(gene.info.table, mean.svm, by = setNames("GENE_entrezID", gene.column.name)) %>% 
      inner_join(., mean.fp, by = setNames("GENE_entrezID", gene.column.name)) %>% 
      clean_names()
    
    saveRDS(final.table, results.file)
    return(final.table)

    all.results[[i]] <- final.table
   
	}
}
  
    
