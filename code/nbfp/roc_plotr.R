roc_plotr <- function(score_tab = NULL, 
                      gene_col = "gene", 
                      pval_col = "logGWASp", 
                      comp_score_col = "logFDR", 
                      comb_score_col = "comb_score", 
                      disgenet_genes = NULL, 
                      conf_cutoffs = c(0.2, 0.3, 0.4, 0.5)) {
  
  score_tab_renamed <- score_tab %>% 
    rename(gene = {{ gene_col }}, log_pval = {{ pval_col }}, comp_score = {{ comp_score_col }}, comb_score = {{ comb_score_col }})
  
  pval_rocs <- list()
  for(i in 1:length(conf_cutoffs)) {
    cutoff_evidence <- disgenet_genes %>% 
      filter(score > conf_cutoffs[i]) %>% 
      dplyr::select(gene_symbol)
    
    score_tab_labelled <- score_tab_renamed %>% 
      mutate(gwas_sig = ifelse(gene %in% cutoff_evidence$gene_symbol, TRUE, FALSE))
    
    pval_rocs[[i]] <- roc(score_tab_labelled$gwas_sig, score_tab_labelled$log_pval, levels = c(TRUE, FALSE))
    
  }
  
  names(pval_rocs) <- conf_cutoffs
  pval_roc_plots <- ggroc(pval_rocs, legacy.axes = TRUE) +
    ggtitle(paste(pval_col, "ROCs at 4 different DisGeNET Confidence Levels")) +
    geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color="grey", linetype="dashed") +
    labs(xlab = "FPR",
         ylab = "TPR",
         color = "Confidence Scores")
  
  comp_score_col_rocs <- list()
  for(i in 1:length(conf_cutoffs)) {
    cutoff_evidence <- disgenet_genes %>% 
      filter(score > conf_cutoffs[i]) %>% 
      dplyr::select(gene_symbol)
    
    score_tab_labelled <- score_tab_renamed %>% 
      mutate(gwas_sig = ifelse(gene %in% cutoff_evidence$gene_symbol, TRUE, FALSE))
    
    comp_score_col_rocs[[i]] <- roc(score_tab_labelled$gwas_sig, score_tab_labelled$comp_score, levels = c(TRUE, FALSE))
    
  }
  
  names(comp_score_col_rocs) <- conf_cutoffs
  comp_score_col_roc_plots <- ggroc(comp_score_col_rocs, legacy.axes = TRUE) +
    ggtitle(paste("SIGNET", comp_score_col, "ROCs at 4 different DisGeNET Confidence Levels")) +
    geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color="grey", linetype="dashed") +
    labs(xlab = "FPR",
         ylab = "TPR",
         color = "Confidence Scores")
  
  
  comb_score_rocs <- list()
  for(i in 1:length(conf_cutoffs)) {
    cutoff_evidence <- disgenet_genes %>% 
      filter(score > conf_cutoffs[i]) %>% 
      dplyr::select(gene_symbol)
    
    score_tab_labelled <- score_tab_renamed %>% 
      mutate(gwas_sig = ifelse(gene %in% cutoff_evidence$gene_symbol, TRUE, FALSE))
    
    comb_score_rocs[[i]] <- roc(score_tab_labelled$gwas_sig, score_tab_labelled$comb_score, levels = c(TRUE, FALSE))
    
  }
  
  names(comb_score_rocs) <- conf_cutoffs
  comb_score_roc_plots <- ggroc(comb_score_rocs, legacy.axes = TRUE) +
    ggtitle(paste(comb_score_col, "ROCs at 4 different DisGeNET Confidence Levels")) +
    geom_segment(aes(x = 1, xend = 0, y = 1, yend = 0), color="grey", linetype="dashed") +
    labs(xlab = "FPR",
         ylab = "TPR",
         color = "Confidence Scores")
  
  comb_score_roc_plots + comp_score_col_roc_plots + pval_roc_plots
    
}