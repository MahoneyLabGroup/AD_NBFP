#Pulls in table made from merge.svm.gene.info.manager() and adds -log10 transformed FPR and Pval columns
#Performs a gene p-value cutoff, default is 10e-3
make.table.for.pareto <- function(results.dir, p_cutoff = 10e-3, normalize = FALSE, do_cutoff = TRUE) {
  
  module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
  module.dir <- module.dir.info$module.dir
  dir.table <- module.dir.info$dir.table
  
  all.results <- vector(mode = "list", length = length(module.dir))
  names(all.results) <- dir.table[,2]
  
  for(i in 1:length(module.dir)){
    
    #Read in final.table from merge.svm.gene.info.manager
    fin.tab.path <- paste0(module.dir[i], "/Candidate.Gene.Results.csv")
    fin.tab <- read_csv(fin.tab.path, col_names = TRUE)
    
    #Add two -log10 transformed columns (P and FPR)
    fin.tab$logFP <- -log10(as.numeric(fin.tab$Mean.FP.Rate))
    fin.tab$logPV <- -log10(as.numeric(fin.tab$P))
    
    fin.tab$BetterThan <- 0 
    #Calculate new ranking for each gene based on log transformed data
    for(i in 1:length(fin.tab$GENE)){
      fin.tab$BetterThan[[i]] <- length(intersect(fin.tab$GENE[fin.tab$logFP < fin.tab$logFP[[i]]],
                                                  fin.tab$GENE[fin.tab$logPV < fin.tab$logPV[[i]]]))
    }
    
    #Subset data to only take significant genes to plot
    if(do_cutoff){
    fin.tab %<>% dplyr::filter(P < p_cutoff)
    }
    
    #if normalize is TRUE perform normalization
    if(normalize) {
      fin.tab.norm <- fin.tab
      fin.tab.norm$logFP <- fin.tab.norm$logFP/max(fin.tab.norm$logFP)
      fin.tab.norm$logPV <- fin.tab.norm$logPV/max(fin.tab.norm$logPV)
      
      return(fin.tab.norm)
      
    }else{
      
      return(fin.tab)
      
    }
  }
}