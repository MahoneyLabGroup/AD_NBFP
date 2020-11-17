merge_svm_results <- function(results_dir = NULL, gene_info_table, gene_col = "gene") {
  
  results <- list.files(results_dir, full.names = TRUE)
  
  for(i in 1:length(results)){
    
    tab <- read_csv(results[[i]], col_names = TRUE) %>% 
      rename(gene = {{ gene_col }}) %>% 
      mutate(across(gene, as.character))
      
    
  }

}