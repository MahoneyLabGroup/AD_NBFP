gene_results_formatr <- function(score_tab = NULL,
                                 gene_id = "gene",
                                 gene_name = "name",
                                 pval_col = "p",
                                 fpr_col = "mean_fp_rate",
                                 make_scaled_vals = FALSE) {
  score_tab_renamed <- score_tab %>% 
    rename(gene_id = {{ gene_id }}, gene_name = {{ gene_name }}, pval = {{ pval_col }}, fpr = {{ fpr_col }}) %>%
    mutate(
      log_p = -log10(pval),
      log_fpr = -log10(fpr) 
    ) %>% 
    rowwise() %>%
    mutate(comb_score = sum(log_p > .$log_p & log_fpr > .$log_fpr)) %>%
    ungroup() %>%
    mutate(comb_score = comb_score/n()) %>% 
    select(gene_id, 
           gene_name, 
           log_p, 
           log_fpr,
           comb_score,
           chromosome_name,
           mean_gene_position,
           start_position,
           end_position)
  
  if(make_scaled_vals){
    score_tab_renamed %<>% 
      mutate(
        log_p_scaled_to_1 = log_p/max(log_p),
        log_fpr_scaled_to_1 = log_fpr/max(log_fpr)
      )
  }
  return(score_tab_renamed)
}