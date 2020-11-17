auc_plotr <- function(score_tab = NULL, 
                      gene_col = "gene", 
                      pval_col = "logGWASp", 
                      comp_score_col = "logFDR", 
                      comb_score_col = "comb_score", 
                      disgenet_genes = NULL, 
                      conf_cutoffs = c(0.2, 0.3, 0.4, 0.5),
                      title = "AUC Comparisons") {
  
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
  pval_auc_table <- tribble(
    ~AUC, ~conf_score, ~col_type,
    as.numeric(pval_rocs[[1]]$auc), conf_cutoffs[[1]], "pval",
    as.numeric(pval_rocs[[2]]$auc), conf_cutoffs[[2]], "pval",
    as.numeric(pval_rocs[[3]]$auc), conf_cutoffs[[3]], "pval",
    as.numeric(pval_rocs[[4]]$auc), conf_cutoffs[[4]], "pval",
  )
  
  
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
  comp_score_col_auc_table <- tribble(
    ~AUC, ~conf_score, ~col_type,
    as.numeric(comp_score_col_rocs[[1]]$auc), conf_cutoffs[[1]], "comp_score",
    as.numeric(comp_score_col_rocs[[2]]$auc), conf_cutoffs[[2]], "comp_score",
    as.numeric(comp_score_col_rocs[[3]]$auc), conf_cutoffs[[3]], "comp_score",
    as.numeric(comp_score_col_rocs[[4]]$auc), conf_cutoffs[[4]], "comp_score",
  )
  
  
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
  comb_score_auc_table <- tribble(
    ~AUC, ~conf_score, ~col_type,
    as.numeric(comb_score_rocs[[1]]$auc), conf_cutoffs[[1]], "comb_score",
    as.numeric(comb_score_rocs[[2]]$auc), conf_cutoffs[[2]], "comb_score",
    as.numeric(comb_score_rocs[[3]]$auc), conf_cutoffs[[3]], "comb_score",
    as.numeric(comb_score_rocs[[4]]$auc), conf_cutoffs[[4]], "comb_score",
  )
  
  
  all_aucs <- bind_rows(pval_auc_table, comp_score_col_auc_table, comb_score_auc_table)
    
  auc_comp_plot <- all_aucs %>% 
    ggplot(aes(x = reorder_within(conf_score, AUC, col_type), y = AUC, fill = conf_score)) +
    geom_col() +
    facet_wrap(~col_type, scales = "free", nrow = 1, ncol = 3) +
    coord_flip() +
    scale_x_reordered() +
    labs(x = "",
         title = title) +
    expand_limits(y = 1) +
    theme_light() +
    theme(legend.position = "none")
  
}