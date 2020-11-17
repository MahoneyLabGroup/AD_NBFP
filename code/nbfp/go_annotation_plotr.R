go_annotation_plotr <- function(score_tab = NULL, 
                                gene_col = "gene_name",
                                gwas_p_scores = "log_p",
                                comparison_scores = "log_fpr"){

    score_tab_renamed <- score_tab %>% 
      rename(gene_name = {{ gene_col }}, p_score = {{ gwas_p_scores }}, comp_scores = {{ comparison_scores }}) %>% 
      rowwise() %>%
      mutate(comb_score = sum(p_score > .$p_score & comp_scores > .$comp_scores)) %>%
      ungroup() %>%
      mutate(comb_score = comb_score/n())

  
  #Ontology for P value
  pval_gsea <- score_tab_renamed %>%  
    arrange(desc(p_score))
  
  pval_gsea_gost <- gprofiler2::gost(pval_gsea$gene_name, 
                                     organism = "hsapiens",
                                     ordered_query = TRUE,
                                     numeric_ns = "ENTREZGENE_ACC",
                                     sources = c("GO:BP", "REAC", "KEGG"))
  
  pval_gsea_res <- pval_gsea_gost$result %>% 
    select(term_id, p_val_gsea_term = term_name, p_val_gsea_p = p_value)
  
  #Ontology for combined score
  comb_score_gsea <- score_tab_renamed %>% 
    arrange(desc(comb_score))
  
  comb_score_gsea_gost <- gprofiler2::gost(comb_score_gsea$gene_name,
                                      organism = "hsapiens",
                                      ordered_query = TRUE,
                                      numeric_ns = "ENTREZGENE_ACC",
                                      sources = c("GO:BP", "REAC", "KEGG"))
  
  comb_score_gsea_res <- comb_score_gsea_gost$result %>% 
    select(term_id, comb_score_gsea_term = term_name, comb_score_gsea_p = p_value)
  
  #Combining tables and plotting
  combined_res <- pval_gsea_res %>% full_join(comb_score_gsea_res, by = "term_id") %>% 
    mutate(
      log_p_p = -log10(p_val_gsea_p),
      log_comb_score_p = -log10(comb_score_gsea_p)
    ) %>% 
    mutate(across(starts_with("log_"), replace_na, 0)) %>% 
    mutate(
      p_val_gsea_term = ifelse(is.na(p_val_gsea_term), comb_score_gsea_term, p_val_gsea_term),
      comb_score_gsea_term = ifelse(is.na(comb_score_gsea_term), p_val_gsea_term, comb_score_gsea_term)
    )
  
  
  plot <- combined_res %>% 
    ggplot(aes(log_p_p, log_comb_score_p)) +
    geom_point() +
    geom_abline(slope = 1, 
                intercept = 0) +
    geom_text_repel(aes(label = ifelse(log_p_p > 0 & log_comb_score_p > 0, as.character(comb_score_gsea_term), "")),
              hjust = -0.2, 
              vjust = -0.5, 
              angle = 15, 
              size = 2) +
    geom_text(aes(label = as.character(p_val_gsea_term)), data = combined_res %>% filter(log_comb_score_p == 0) %>% arrange(desc(log_p_p)) %>%  slice_head(n = 5),
              hjust = 0.2, 
              vjust = 0.5, 
              angle = 345, 
              size = 2) +
    geom_text(aes(label = as.character(comb_score_gsea_term)), data = combined_res %>% filter(log_p_p == 0) %>% arrange(desc(log_comb_score_p)) %>% slice_head(n = 5),
              hjust = 0.2, 
              vjust = 0.5, 
              angle = 0, 
              size = 2) +
    #expand_limits() +
    labs(x = "-log10 of Pvalue for GWAS P",
         y = "-log10 of Pvalue for Combined Score")
  
  return(plot)
}